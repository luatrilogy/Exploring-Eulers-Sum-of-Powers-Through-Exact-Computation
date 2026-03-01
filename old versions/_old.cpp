#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <boost/multiprecision/cpp_int.hpp>

#include <algorithm>
#include <cstdint>
#include <stdexcept>
#include <unordered_set>
#include <vector>

namespace py = pybind11;
using boost::multiprecision::cpp_int;
// Portable fixed-width 128-bit unsigned type (works on MSVC too)
using u128 = boost::multiprecision::uint128_t;

// -------------------- math helpers --------------------

static inline uint32_t powmod_u32(uint32_t a, int k, uint32_t mod) {
    uint64_t res = 1;
    uint64_t base = (uint64_t)(a % mod);
    int e = k;
    while (e > 0) {
        if (e & 1) res = (res * base) % mod;
        base = (base * base) % mod;
        e >>= 1;
    }
    return (uint32_t)res;
}

static inline cpp_int pow_cpp_u32(uint32_t a, int k) {
    cpp_int res = 1;
    cpp_int base = a;
    int e = k;
    while (e > 0) {
        if (e & 1) res *= base;
        base *= base;
        e >>= 1;
    }
    return res;
}

static inline u128 pow_u128_fast(uint32_t a, int k) {
    u128 res = 1;
    u128 base = (u128)a;
    int e = k;
    while (e > 0) {
        if (e & 1) res *= base;
        e >>= 1;
        if (e) base *= base;
    }
    return res;
}

static inline uint32_t gcd_u32(uint32_t a, uint32_t b) {
    while (b) { uint32_t t = a % b; a = b; b = t; }
    return a;
}

static inline uint64_t comb_u64(uint64_t n, uint64_t r) {
    if (r > n) return 0;
    if (r > n - r) r = n - r;
    uint64_t num = 1, den = 1;
    for (uint64_t i = 1; i <= r; ++i) {
        // exact division at end is safe for these small r (3 or 4)
        num *= (n - r + i);
        den *= i;
    }
    return num / den;
}

static inline void progress_call(const py::object& progress_cb, uint64_t done, uint64_t total) {
    if (progress_cb.is_none()) return;
    py::gil_scoped_acquire gil;
    progress_cb((uint64_t)done, (uint64_t)total);
}

// -------------------- packing + witnesses --------------------

struct W3u16 { uint16_t a,b,c; };
struct W3u32 { uint32_t a,b,c; };
struct W4u16 { uint16_t a,b,c,d; };
struct W4u32 { uint32_t a,b,c,d; };

// pack two residues (< 65536) into 32-bit key
static inline uint32_t pack_key16(uint16_t r1, uint16_t r2) {
    return (uint32_t(r1) << 16) | uint32_t(r2);
}

// Compute gcd of all terms AND b
static inline uint32_t gcd_terms_and_b(const std::vector<uint32_t>& terms, uint32_t b) {
    uint32_t g = b;
    for (uint32_t x : terms) g = gcd_u32(g, x);
    return g;
}

// very small signature for de-dup after primitive grouping
static inline uint64_t signature_terms_b(const std::vector<uint32_t>& terms, uint32_t b) {
    uint64_t h = 1469598103934665603ull;
    for (uint32_t v : terms) {
        h ^= (uint64_t)v + 0x9e3779b97f4a7c15ull + (h<<6) + (h>>2);
    }
    h ^= (uint64_t)b + 0x9e3779b97f4a7c15ull + (h<<6) + (h>>2);
    return h;
}

template <typename EntryT>
static inline void sort_bucket_by_key(std::vector<EntryT>& v) {
    std::sort(v.begin(), v.end(), [](const EntryT& x, const EntryT& y){ return x.key < y.key; });
}

// Decide whether all needed exact values fit into u128 for fast cutoff/verify.
// We test max_val^k and max_b^k using cpp_int and compare to 2^128.
static inline bool can_use_u128(uint32_t max_val, uint32_t max_b, int k) {
    cpp_int limit = cpp_int(1);
    limit <<= 128;
    cpp_int a = pow_cpp_u32(max_val, k);
    cpp_int b = pow_cpp_u32(max_b, k);
    // also need up to 8 terms sum; if each fits, sum fits if 8*a < 2^128
    return (a < limit) && (b < limit) && (a * 8 < limit);
}

// -------------------- m=8 split 4+4 --------------------

struct Entry4u16 { uint32_t key; W4u16 w; };
struct Entry4u32 { uint32_t key; W4u32 w; };

static py::list search_m8_blocked_sieved_impl(
    int k,
    uint32_t min_b,
    uint32_t max_b,
    uint32_t max_val,
    uint32_t mod_M,
    uint32_t p1,
    uint32_t p2,
    uint32_t block_b,
    uint32_t max_solutions,
    uint64_t max_entries_cap,
    bool primitive_grouping,
    py::object progress_cb
) {
    if (k <= 0) throw std::invalid_argument("k must be positive");
    if (min_b < 2) min_b = 2;
    if (max_b < min_b) return py::list();

    if (max_val == 0) max_val = max_b; // match Python default: terms up to max_b
    if (max_val < 2) max_val = 2;

    // For packed mod arithmetic we require power-of-two mod_M = 2^16 to use uint16 wrap.
    if (mod_M != 65536u) throw std::invalid_argument("mod_M must be 65536 (2^16) in this backend");
    if (p1 >= 65536u || p2 >= 65536u) throw std::invalid_argument("p1 and p2 must be < 65536");

    const uint64_t total_b = (uint64_t)max_b - (uint64_t)min_b + 1;

    // #nondecreasing 4-tuples from 1..max_val is C(max_val+3,4)
    const uint64_t est_entries = comb_u64((uint64_t)max_val + 3, 4);
    if (max_entries_cap != 0 && est_entries > max_entries_cap) {
        throw std::runtime_error("LHS table would exceed max_entries_cap. Reduce max_val or raise the cap.");
    }

    // Residues for 0..max_val
    std::vector<uint16_t> powM(max_val + 1);
    std::vector<uint16_t> powP1(max_val + 1);
    std::vector<uint16_t> powP2(max_val + 1);
    for (uint32_t i = 0; i <= max_val; ++i) {
        powM[i]  = (uint16_t)powmod_u32(i, k, mod_M);
        powP1[i] = (uint16_t)powmod_u32(i, k, p1);
        powP2[i] = (uint16_t)powmod_u32(i, k, p2);
    }

    const bool use_u16 = (max_val <= 65535u);
    const bool use_u128 = can_use_u128(max_val, max_b, k);

    std::vector<u128> powExact128;
    if (use_u128) {
        powExact128.resize(max_val + 1);
        for (uint32_t i = 0; i <= max_val; ++i) {
            powExact128[i] = pow_u128_fast(i, k);
        }
    }

    std::vector<u128> powB128;
    if (use_u128) {
        powB128.resize(max_b + 1);
        for (uint32_t i = 0; i <= max_b; ++i) powB128[i] = pow_u128_fast(i, k);
    }

    auto verify_u128 = [&](const std::vector<uint32_t>& terms, uint32_t B) -> bool {
        u128 s = 0;
        for (uint32_t x : terms) s += powExact128[x];
        return s == powB128[B];
    };

    auto verify_cpp = [&](const std::vector<uint32_t>& terms, uint32_t B) -> bool {
        cpp_int s = 0;
        for (uint32_t x : terms) s += pow_cpp_u32(x, k);
        return s == pow_cpp_u32(B, k);
    };

    py::list out;
    std::unordered_set<uint64_t> seen;
    if (primitive_grouping) seen.reserve((size_t)max_solutions * 2);

    if (use_u16) {
        // LHS bucketed by mod_M; within each bucket sorted by (p1,p2) key.
        std::vector<std::vector<Entry4u16>> buckets(mod_M);
        const uint64_t avg = est_entries / mod_M + 8;
        for (auto& v : buckets) v.reserve((size_t)avg);

        uint64_t total_entries = 0;
        for (uint32_t a = 1; a <= max_val; ++a) {
            const uint16_t aM = powM[a], a1 = powP1[a], a2 = powP2[a];
            for (uint32_t b = a; b <= max_val; ++b) {
                const uint16_t abM = (uint16_t)(aM + powM[b]);
                uint32_t ab1 = (uint32_t)a1 + (uint32_t)powP1[b]; if (ab1 >= p1) ab1 -= p1;
                uint32_t ab2 = (uint32_t)a2 + (uint32_t)powP2[b]; if (ab2 >= p2) ab2 -= p2;
                for (uint32_t c = b; c <= max_val; ++c) {
                    const uint16_t abcM = (uint16_t)(abM + powM[c]);
                    uint32_t abc1 = ab1 + (uint32_t)powP1[c]; if (abc1 >= p1) abc1 -= p1;
                    uint32_t abc2 = ab2 + (uint32_t)powP2[c]; if (abc2 >= p2) abc2 -= p2;
                    for (uint32_t d = c; d <= max_val; ++d) {
                        const uint16_t sM = (uint16_t)(abcM + powM[d]);
                        uint32_t s1 = abc1 + (uint32_t)powP1[d]; if (s1 >= p1) s1 -= p1;
                        uint32_t s2 = abc2 + (uint32_t)powP2[d]; if (s2 >= p2) s2 -= p2;

                        buckets[sM].push_back(Entry4u16{pack_key16((uint16_t)s1, (uint16_t)s2),
                                                       W4u16{(uint16_t)a,(uint16_t)b,(uint16_t)c,(uint16_t)d}});
                        total_entries++;
                        if (max_entries_cap != 0 && total_entries > max_entries_cap) {
                            throw std::runtime_error("Exceeded max_entries_cap while building LHS table");
                        }
                    }
                }
            }
        }
        for (auto& v : buckets) sort_bucket_by_key(v);

        // Scan B
        uint64_t done = 0;
        for (uint32_t bStart = min_b; bStart <= max_b; ) {
            uint32_t bEnd = bStart + block_b - 1;
            if (bEnd > max_b) bEnd = max_b;

            for (uint32_t B = bStart; B <= bEnd; ++B) {
                done++;
                progress_call(progress_cb, done, total_b);

                const uint16_t tM  = (uint16_t)powmod_u32(B, k, mod_M);
                const uint16_t tP1 = (uint16_t)powmod_u32(B, k, p1);
                const uint16_t tP2 = (uint16_t)powmod_u32(B, k, p2);

                const uint32_t rhs_limit = (B < max_val ? B : max_val);

                u128 target128 = 0;
                if (use_u128) target128 = powB128[B];

                for (uint32_t a = 1; a <= rhs_limit; ++a) {
                    const uint16_t aM = powM[a], a1 = powP1[a], a2 = powP2[a];
                    for (uint32_t bb = a; bb <= rhs_limit; ++bb) {
                        const uint16_t abM = (uint16_t)(aM + powM[bb]);
                        uint32_t ab1 = (uint32_t)a1 + (uint32_t)powP1[bb]; if (ab1 >= p1) ab1 -= p1;
                        uint32_t ab2 = (uint32_t)a2 + (uint32_t)powP2[bb]; if (ab2 >= p2) ab2 -= p2;
                        for (uint32_t c = bb; c <= rhs_limit; ++c) {
                            const uint16_t abcM = (uint16_t)(abM + powM[c]);
                            uint32_t abc1 = ab1 + (uint32_t)powP1[c]; if (abc1 >= p1) abc1 -= p1;
                            uint32_t abc2 = ab2 + (uint32_t)powP2[c]; if (abc2 >= p2) abc2 -= p2;

                            u128 abc128 = 0;
                            if (use_u128) abc128 = powExact128[a] + powExact128[bb] + powExact128[c];

                            for (uint32_t d = c; d <= rhs_limit; ++d) {
                                if (use_u128) {
                                    u128 rhs4 = abc128 + powExact128[d];
                                    if (rhs4 > target128) break; // safe since d increases
                                }

                                const uint16_t rM = (uint16_t)(abcM + powM[d]);
                                uint32_t r1 = abc1 + (uint32_t)powP1[d]; if (r1 >= p1) r1 -= p1;
                                uint32_t r2 = abc2 + (uint32_t)powP2[d]; if (r2 >= p2) r2 -= p2;

                                const uint16_t needM  = (uint16_t)(tM - rM);
                                const uint16_t needP1 = (uint16_t)((tP1 + p1 - (uint16_t)r1) % p1);
                                const uint16_t needP2 = (uint16_t)((tP2 + p2 - (uint16_t)r2) % p2);
                                const uint32_t needKey = pack_key16(needP1, needP2);

                                const auto& vec = buckets[needM];
                                if (vec.empty()) continue;

                                auto it = std::lower_bound(vec.begin(), vec.end(), needKey,
                                    [](const Entry4u16& e, uint32_t key){ return e.key < key; });
                                for (; it != vec.end() && it->key == needKey; ++it) {
                                    std::vector<uint32_t> terms = {
                                        (uint32_t)it->w.a,(uint32_t)it->w.b,(uint32_t)it->w.c,(uint32_t)it->w.d,
                                        a,bb,c,d
                                    };
                                    std::sort(terms.begin(), terms.end());

                                    bool ok = use_u128 ? verify_u128(terms, B) : verify_cpp(terms, B);
                                    if (!ok) continue;

                                    uint32_t outB = B;
                                    if (primitive_grouping) {
                                        uint32_t g = gcd_terms_and_b(terms, B);
                                        if (g > 1) {
                                            for (auto& v : terms) v /= g;
                                            outB /= g;
                                        }
                                        uint64_t sig = signature_terms_b(terms, outB);
                                        if (seen.find(sig) != seen.end()) continue;
                                        seen.insert(sig);
                                    }

                                    py::tuple lhs8(8);
                                    for (int i = 0; i < 8; ++i) lhs8[i] = terms[i];
                                    out.append(py::make_tuple(lhs8, outB));
                                    if ((uint32_t)out.size() >= max_solutions) return out;
                                }
                            }
                        }
                    }
                }
            }
            bStart = bEnd + 1;
        }
        return out;
    }

    // Fallback (max_val > 65535): use u32 witnesses
    {
        std::vector<std::vector<Entry4u32>> buckets(mod_M);
        const uint64_t avg = est_entries / mod_M + 8;
        for (auto& v : buckets) v.reserve((size_t)avg);

        uint64_t total_entries = 0;
        for (uint32_t a = 1; a <= max_val; ++a) {
            const uint16_t aM = powM[a], a1 = powP1[a], a2 = powP2[a];
            for (uint32_t b = a; b <= max_val; ++b) {
                const uint16_t abM = (uint16_t)(aM + powM[b]);
                uint32_t ab1 = (uint32_t)a1 + (uint32_t)powP1[b]; if (ab1 >= p1) ab1 -= p1;
                uint32_t ab2 = (uint32_t)a2 + (uint32_t)powP2[b]; if (ab2 >= p2) ab2 -= p2;
                for (uint32_t c = b; c <= max_val; ++c) {
                    const uint16_t abcM = (uint16_t)(abM + powM[c]);
                    uint32_t abc1 = ab1 + (uint32_t)powP1[c]; if (abc1 >= p1) abc1 -= p1;
                    uint32_t abc2 = ab2 + (uint32_t)powP2[c]; if (abc2 >= p2) abc2 -= p2;
                    for (uint32_t d = c; d <= max_val; ++d) {
                        const uint16_t sM = (uint16_t)(abcM + powM[d]);
                        uint32_t s1 = abc1 + (uint32_t)powP1[d]; if (s1 >= p1) s1 -= p1;
                        uint32_t s2 = abc2 + (uint32_t)powP2[d]; if (s2 >= p2) s2 -= p2;

                        buckets[sM].push_back(Entry4u32{pack_key16((uint16_t)s1, (uint16_t)s2), W4u32{a,b,c,d}});
                        total_entries++;
                        if (max_entries_cap != 0 && total_entries > max_entries_cap) {
                            throw std::runtime_error("Exceeded max_entries_cap while building LHS table");
                        }
                    }
                }
            }
        }
        for (auto& v : buckets) sort_bucket_by_key(v);

        uint64_t done = 0;
        for (uint32_t bStart = min_b; bStart <= max_b; ) {
            uint32_t bEnd = bStart + block_b - 1;
            if (bEnd > max_b) bEnd = max_b;

            for (uint32_t B = bStart; B <= bEnd; ++B) {
                done++;
                progress_call(progress_cb, done, total_b);

                const uint16_t tM  = (uint16_t)powmod_u32(B, k, mod_M);
                const uint16_t tP1 = (uint16_t)powmod_u32(B, k, p1);
                const uint16_t tP2 = (uint16_t)powmod_u32(B, k, p2);

                const uint32_t rhs_limit = (B < max_val ? B : max_val);

                for (uint32_t a = 1; a <= rhs_limit; ++a) {
                    const uint16_t aM = powM[a], a1 = powP1[a], a2 = powP2[a];
                    for (uint32_t bb = a; bb <= rhs_limit; ++bb) {
                        const uint16_t abM = (uint16_t)(aM + powM[bb]);
                        uint32_t ab1 = (uint32_t)a1 + (uint32_t)powP1[bb]; if (ab1 >= p1) ab1 -= p1;
                        uint32_t ab2 = (uint32_t)a2 + (uint32_t)powP2[bb]; if (ab2 >= p2) ab2 -= p2;
                        for (uint32_t c = bb; c <= rhs_limit; ++c) {
                            const uint16_t abcM = (uint16_t)(abM + powM[c]);
                            uint32_t abc1 = ab1 + (uint32_t)powP1[c]; if (abc1 >= p1) abc1 -= p1;
                            uint32_t abc2 = ab2 + (uint32_t)powP2[c]; if (abc2 >= p2) abc2 -= p2;
                            for (uint32_t d = c; d <= rhs_limit; ++d) {
                                const uint16_t rM = (uint16_t)(abcM + powM[d]);
                                uint32_t r1 = abc1 + (uint32_t)powP1[d]; if (r1 >= p1) r1 -= p1;
                                uint32_t r2 = abc2 + (uint32_t)powP2[d]; if (r2 >= p2) r2 -= p2;

                                const uint16_t needM  = (uint16_t)(tM - rM);
                                const uint16_t needP1 = (uint16_t)((tP1 + p1 - (uint16_t)r1) % p1);
                                const uint16_t needP2 = (uint16_t)((tP2 + p2 - (uint16_t)r2) % p2);
                                const uint32_t needKey = pack_key16(needP1, needP2);

                                const auto& vec = buckets[needM];
                                if (vec.empty()) continue;

                                auto it = std::lower_bound(vec.begin(), vec.end(), needKey,
                                    [](const Entry4u32& e, uint32_t key){ return e.key < key; });
                                for (; it != vec.end() && it->key == needKey; ++it) {
                                    std::vector<uint32_t> terms = {it->w.a,it->w.b,it->w.c,it->w.d, a,bb,c,d};
                                    std::sort(terms.begin(), terms.end());
                                    if (!verify_cpp(terms, B)) continue;

                                    uint32_t outB = B;
                                    if (primitive_grouping) {
                                        uint32_t g = gcd_terms_and_b(terms, B);
                                        if (g > 1) {
                                            for (auto& v : terms) v /= g;
                                            outB /= g;
                                        }
                                        uint64_t sig = signature_terms_b(terms, outB);
                                        if (seen.find(sig) != seen.end()) continue;
                                        seen.insert(sig);
                                    }

                                    py::tuple lhs8(8);
                                    for (int i = 0; i < 8; ++i) lhs8[i] = terms[i];
                                    out.append(py::make_tuple(lhs8, outB));
                                    if ((uint32_t)out.size() >= max_solutions) return out;
                                }
                            }
                        }
                    }
                }
            }
            bStart = bEnd + 1;
        }
        return out;
    }
}

// -------------------- m=6 split 3+3 --------------------

struct Entry3u16 { uint32_t key; W3u16 w; };
struct Entry3u32 { uint32_t key; W3u32 w; };

static py::list search_m6_3x3_impl(
    int k,
    uint32_t min_b,
    uint32_t max_b,
    uint32_t max_val,
    uint32_t mod_M,
    uint32_t p1,
    uint32_t p2,
    uint32_t block_b,
    uint32_t max_solutions,
    uint64_t max_entries_cap,
    bool primitive_grouping,
    py::object progress_cb
) {
    if (k <= 0) throw std::invalid_argument("k must be positive");
    if (min_b < 2) min_b = 2;
    if (max_b < min_b) return py::list();

    if (max_val == 0) max_val = max_b;
    if (max_val < 2) max_val = 2;

    if (mod_M != 65536u) throw std::invalid_argument("mod_M must be 65536 (2^16) in this backend");
    if (p1 >= 65536u || p2 >= 65536u) throw std::invalid_argument("p1 and p2 must be < 65536");

    const uint64_t total_b = (uint64_t)max_b - (uint64_t)min_b + 1;

    // #nondecreasing 3-tuples from 1..max_val is C(max_val+2,3)
    const uint64_t est_entries = comb_u64((uint64_t)max_val + 2, 3);
    if (max_entries_cap != 0 && est_entries > max_entries_cap) {
        throw std::runtime_error("LHS table would exceed max_entries_cap. Reduce max_val or raise the cap.");
    }

    std::vector<uint16_t> powM(max_val + 1);
    std::vector<uint16_t> powP1(max_val + 1);
    std::vector<uint16_t> powP2(max_val + 1);
    for (uint32_t i = 0; i <= max_val; ++i) {
        powM[i]  = (uint16_t)powmod_u32(i, k, mod_M);
        powP1[i] = (uint16_t)powmod_u32(i, k, p1);
        powP2[i] = (uint16_t)powmod_u32(i, k, p2);
    }

    const bool use_u16 = (max_val <= 65535u);
    const bool use_u128 = can_use_u128(max_val, max_b, k);

    std::vector<u128> powExact128;
    if (use_u128) {
        powExact128.resize(max_val + 1);
        for (uint32_t i = 0; i <= max_val; ++i) {
            powExact128[i] = pow_u128_fast(i, k);
        }
    }

    std::vector<u128> powB128;
    if (use_u128) {
        powB128.resize(max_b + 1);
        for (uint32_t i = 0; i <= max_b; ++i) powB128[i] = pow_u128_fast(i, k);
    }

    auto verify_u128 = [&](const std::vector<uint32_t>& terms, uint32_t B) -> bool {
        u128 s = 0;
        for (uint32_t x : terms) s += powExact128[x];
        return s == powB128[B];
    };
    auto verify_cpp = [&](const std::vector<uint32_t>& terms, uint32_t B) -> bool {
        cpp_int s = 0;
        for (uint32_t x : terms) s += pow_cpp_u32(x, k);
        return s == pow_cpp_u32(B, k);
    };

    py::list out;
    std::unordered_set<uint64_t> seen;
    if (primitive_grouping) seen.reserve((size_t)max_solutions * 2);

    if (use_u16) {
        std::vector<std::vector<Entry3u16>> buckets(mod_M);
        const uint64_t avg = est_entries / mod_M + 8;
        for (auto& v : buckets) v.reserve((size_t)avg);

        uint64_t total_entries = 0;
        for (uint32_t a = 1; a <= max_val; ++a) {
            const uint16_t aM = powM[a], a1 = powP1[a], a2 = powP2[a];
            for (uint32_t b = a; b <= max_val; ++b) {
                const uint16_t abM = (uint16_t)(aM + powM[b]);
                uint32_t ab1 = (uint32_t)a1 + (uint32_t)powP1[b]; if (ab1 >= p1) ab1 -= p1;
                uint32_t ab2 = (uint32_t)a2 + (uint32_t)powP2[b]; if (ab2 >= p2) ab2 -= p2;
                for (uint32_t c = b; c <= max_val; ++c) {
                    const uint16_t sM = (uint16_t)(abM + powM[c]);
                    uint32_t s1 = ab1 + (uint32_t)powP1[c]; if (s1 >= p1) s1 -= p1;
                    uint32_t s2 = ab2 + (uint32_t)powP2[c]; if (s2 >= p2) s2 -= p2;

                    buckets[sM].push_back(Entry3u16{pack_key16((uint16_t)s1,(uint16_t)s2),
                                                   W3u16{(uint16_t)a,(uint16_t)b,(uint16_t)c}});
                    total_entries++;
                    if (max_entries_cap != 0 && total_entries > max_entries_cap) {
                        throw std::runtime_error("Exceeded max_entries_cap while building LHS table");
                    }
                }
            }
        }
        for (auto& v : buckets) sort_bucket_by_key(v);

        uint64_t done = 0;
        for (uint32_t bStart = min_b; bStart <= max_b; ) {
            uint32_t bEnd = bStart + block_b - 1;
            if (bEnd > max_b) bEnd = max_b;

            for (uint32_t B = bStart; B <= bEnd; ++B) {
                done++;
                progress_call(progress_cb, done, total_b);

                const uint16_t tM  = (uint16_t)powmod_u32(B, k, mod_M);
                const uint16_t tP1 = (uint16_t)powmod_u32(B, k, p1);
                const uint16_t tP2 = (uint16_t)powmod_u32(B, k, p2);

                const uint32_t rhs_limit = (B < max_val ? B : max_val);

                u128 target128 = 0;
                if (use_u128) target128 = powB128[B];

                for (uint32_t d = 1; d <= rhs_limit; ++d) {
                    const uint16_t dM = powM[d], d1 = powP1[d], d2 = powP2[d];
                    for (uint32_t e = d; e <= rhs_limit; ++e) {
                        const uint16_t deM = (uint16_t)(dM + powM[e]);
                        uint32_t de1 = (uint32_t)d1 + (uint32_t)powP1[e]; if (de1 >= p1) de1 -= p1;
                        uint32_t de2 = (uint32_t)d2 + (uint32_t)powP2[e]; if (de2 >= p2) de2 -= p2;

                        u128 de128 = 0;
                        if (use_u128) de128 = powExact128[d] + powExact128[e];

                        for (uint32_t f = e; f <= rhs_limit; ++f) {
                            if (use_u128) {
                                u128 rhs3 = de128 + powExact128[f];
                                if (rhs3 > target128) break;
                            }

                            const uint16_t rM = (uint16_t)(deM + powM[f]);
                            uint32_t r1 = de1 + (uint32_t)powP1[f]; if (r1 >= p1) r1 -= p1;
                            uint32_t r2 = de2 + (uint32_t)powP2[f]; if (r2 >= p2) r2 -= p2;

                            const uint16_t needM  = (uint16_t)(tM - rM);
                            const uint16_t needP1 = (uint16_t)((tP1 + p1 - (uint16_t)r1) % p1);
                            const uint16_t needP2 = (uint16_t)((tP2 + p2 - (uint16_t)r2) % p2);
                            const uint32_t needKey = pack_key16(needP1, needP2);

                            const auto& vec = buckets[needM];
                            if (vec.empty()) continue;

                            auto it = std::lower_bound(vec.begin(), vec.end(), needKey,
                                [](const Entry3u16& e, uint32_t key){ return e.key < key; });
                            for (; it != vec.end() && it->key == needKey; ++it) {
                                std::vector<uint32_t> terms = {
                                    (uint32_t)it->w.a,(uint32_t)it->w.b,(uint32_t)it->w.c,
                                    d,e,f
                                };
                                std::sort(terms.begin(), terms.end());

                                bool ok = use_u128 ? verify_u128(terms, B) : verify_cpp(terms, B);
                                if (!ok) continue;

                                uint32_t outB = B;
                                if (primitive_grouping) {
                                    uint32_t g = gcd_terms_and_b(terms, B);
                                    if (g > 1) {
                                        for (auto& v : terms) v /= g;
                                        outB /= g;
                                    }
                                    uint64_t sig = signature_terms_b(terms, outB);
                                    if (seen.find(sig) != seen.end()) continue;
                                    seen.insert(sig);
                                }

                                py::tuple lhs6(6);
                                for (int i = 0; i < 6; ++i) lhs6[i] = terms[i];
                                out.append(py::make_tuple(lhs6, outB));
                                if ((uint32_t)out.size() >= max_solutions) return out;
                            }
                        }
                    }
                }
            }
            bStart = bEnd + 1;
        }
        return out;
    }

    // Fallback u32 witnesses
    {
        std::vector<std::vector<Entry3u32>> buckets(mod_M);
        const uint64_t avg = est_entries / mod_M + 8;
        for (auto& v : buckets) v.reserve((size_t)avg);

        uint64_t total_entries = 0;
        for (uint32_t a = 1; a <= max_val; ++a) {
            const uint16_t aM = powM[a], a1 = powP1[a], a2 = powP2[a];
            for (uint32_t b = a; b <= max_val; ++b) {
                const uint16_t abM = (uint16_t)(aM + powM[b]);
                uint32_t ab1 = (uint32_t)a1 + (uint32_t)powP1[b]; if (ab1 >= p1) ab1 -= p1;
                uint32_t ab2 = (uint32_t)a2 + (uint32_t)powP2[b]; if (ab2 >= p2) ab2 -= p2;
                for (uint32_t c = b; c <= max_val; ++c) {
                    const uint16_t sM = (uint16_t)(abM + powM[c]);
                    uint32_t s1 = ab1 + (uint32_t)powP1[c]; if (s1 >= p1) s1 -= p1;
                    uint32_t s2 = ab2 + (uint32_t)powP2[c]; if (s2 >= p2) s2 -= p2;

                    buckets[sM].push_back(Entry3u32{pack_key16((uint16_t)s1,(uint16_t)s2), W3u32{a,b,c}});
                    total_entries++;
                    if (max_entries_cap != 0 && total_entries > max_entries_cap) {
                        throw std::runtime_error("Exceeded max_entries_cap while building LHS table");
                    }
                }
            }
        }
        for (auto& v : buckets) sort_bucket_by_key(v);

        uint64_t done = 0;
        for (uint32_t bStart = min_b; bStart <= max_b; ) {
            uint32_t bEnd = bStart + block_b - 1;
            if (bEnd > max_b) bEnd = max_b;

            for (uint32_t B = bStart; B <= bEnd; ++B) {
                done++;
                progress_call(progress_cb, done, total_b);

                const uint16_t tM  = (uint16_t)powmod_u32(B, k, mod_M);
                const uint16_t tP1 = (uint16_t)powmod_u32(B, k, p1);
                const uint16_t tP2 = (uint16_t)powmod_u32(B, k, p2);

                const uint32_t rhs_limit = (B < max_val ? B : max_val);

                for (uint32_t d = 1; d <= rhs_limit; ++d) {
                    const uint16_t dM = powM[d], d1 = powP1[d], d2 = powP2[d];
                    for (uint32_t e = d; e <= rhs_limit; ++e) {
                        const uint16_t deM = (uint16_t)(dM + powM[e]);
                        uint32_t de1 = (uint32_t)d1 + (uint32_t)powP1[e]; if (de1 >= p1) de1 -= p1;
                        uint32_t de2 = (uint32_t)d2 + (uint32_t)powP2[e]; if (de2 >= p2) de2 -= p2;
                        for (uint32_t f = e; f <= rhs_limit; ++f) {
                            const uint16_t rM = (uint16_t)(deM + powM[f]);
                            uint32_t r1 = de1 + (uint32_t)powP1[f]; if (r1 >= p1) r1 -= p1;
                            uint32_t r2 = de2 + (uint32_t)powP2[f]; if (r2 >= p2) r2 -= p2;

                            const uint16_t needM  = (uint16_t)(tM - rM);
                            const uint16_t needP1 = (uint16_t)((tP1 + p1 - (uint16_t)r1) % p1);
                            const uint16_t needP2 = (uint16_t)((tP2 + p2 - (uint16_t)r2) % p2);
                            const uint32_t needKey = pack_key16(needP1, needP2);

                            const auto& vec = buckets[needM];
                            if (vec.empty()) continue;

                            auto it = std::lower_bound(vec.begin(), vec.end(), needKey,
                                [](const Entry3u32& e, uint32_t key){ return e.key < key; });
                            for (; it != vec.end() && it->key == needKey; ++it) {
                                std::vector<uint32_t> terms = {it->w.a,it->w.b,it->w.c, d,e,f};
                                std::sort(terms.begin(), terms.end());
                                if (!verify_cpp(terms, B)) continue;

                                uint32_t outB = B;
                                if (primitive_grouping) {
                                    uint32_t g = gcd_terms_and_b(terms, B);
                                    if (g > 1) {
                                        for (auto& v : terms) v /= g;
                                        outB /= g;
                                    }
                                    uint64_t sig = signature_terms_b(terms, outB);
                                    if (seen.find(sig) != seen.end()) continue;
                                    seen.insert(sig);
                                }

                                py::tuple lhs6(6);
                                for (int i = 0; i < 6; ++i) lhs6[i] = terms[i];
                                out.append(py::make_tuple(lhs6, outB));
                                if ((uint32_t)out.size() >= max_solutions) return out;
                            }
                        }
                    }
                }
            }
            bStart = bEnd + 1;
        }
        return out;
    }
}

// -------------------- pybind11 module --------------------

PYBIND11_MODULE(euler_ext, m) {
    m.doc() = "Euler sum-of-powers MITM core (exact C++ backend)";

    m.def(
        "search_m8_blocked_sieved",
        &search_m8_blocked_sieved_impl,
        py::arg("k"),
        py::arg("min_b"),
        py::arg("max_b"),
        py::arg("max_val") = 0,
        py::arg("mod_M") = 65536,
        py::arg("p1") = 65521,
        py::arg("p2") = 65519,
        py::arg("block_b") = 64,
        py::arg("max_solutions") = 50,
        py::arg("max_entries_cap") = (uint64_t)0, // 0 = no cap; set to guard memory
        py::arg("primitive_grouping") = false,
        py::arg("progress_cb") = py::none()
    );

    m.def(
        "search_m6_3x3",
        &search_m6_3x3_impl,
        py::arg("k"),
        py::arg("min_b"),
        py::arg("max_b"),
        py::arg("max_val") = 0,
        py::arg("mod_M") = 65536,
        py::arg("p1") = 65521,
        py::arg("p2") = 65519,
        py::arg("block_b") = 128,
        py::arg("max_solutions") = 50,
        py::arg("max_entries_cap") = (uint64_t)0,
        py::arg("primitive_grouping") = false,
        py::arg("progress_cb") = py::none()
    );
}