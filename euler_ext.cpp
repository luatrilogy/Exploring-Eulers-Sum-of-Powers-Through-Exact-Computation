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
static inline bool can_use_u128(uint32_t max_val, uint32_t max_b, int k, int max_terms = 8) {
    cpp_int limit = cpp_int(1);
    limit <<= 128;
    cpp_int a = pow_cpp_u32(max_val, k);
    cpp_int b = pow_cpp_u32(max_b, k);
    // Need term^k and b^k to fit, and also (max_terms * term^k) to fit for fast exact checks.
    return (a < limit) && (b < limit) && (a * max_terms < limit);
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
    const bool use_u128 = can_use_u128(max_val, max_b, k, 8);

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
    const bool use_u128 = can_use_u128(max_val, max_b, k, 6);

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


// -------------------- DFS + 2-sum backend (Redesign B) --------------------
//
// Goal: avoid O(n^4) RAM for large max_val by using:
//   * a residue-sieved 2-sum table (O(n^2) memory)
//   * depth-first branch-and-bound for the remaining (m-2) terms
//   * exact verification (u128 when safe, otherwise cpp_int)
//
// This is not "complete" for huge parameter ranges, but it can search domains where 4+4 is impossible.

struct PairEntry {
    uint64_t key;     // packed residues: [modM:16][p1:16][p2:16][p3:16]
    uint16_t a;       // a <= b
    uint16_t b;
    uint16_t maxab;   // == b
};

static inline uint64_t pack_key_4(uint16_t modM, uint16_t r1, uint16_t r2, uint16_t r3) {
    return (uint64_t(modM) << 48) | (uint64_t(r1) << 32) | (uint64_t(r2) << 16) | uint64_t(r3);
}

static inline uint16_t sub_mod_u16(uint16_t t, uint16_t x, uint32_t p) {
    // returns (t - x) mod p, assuming p < 65536
    uint32_t tt = (uint32_t)t;
    uint32_t xx = (uint32_t)x;
    uint32_t res = (tt + p - xx) % p;
    return (uint16_t)res;
}

static inline uint32_t floor_kth_root_u128(u128 val, int k, uint32_t hi) {
    uint32_t lo = 0, r = hi;
    while (lo < r) {
        uint32_t mid = lo + (r - lo + 1) / 2;
        u128 p = pow_u128_fast(mid, k);
        if (p <= val) lo = mid;
        else r = mid - 1;
    }
    return lo;
}

static inline uint32_t ceil_kth_root_u128(u128 val, int k, uint32_t hi) {
    uint32_t x = floor_kth_root_u128(val, k, hi);
    if (pow_u128_fast(x, k) < val) x++;
    return x;
}

static inline uint32_t floor_kth_root_cpp(const cpp_int& val, int k, uint32_t hi) {
    uint32_t lo = 0, r = hi;
    while (lo < r) {
        uint32_t mid = lo + (r - lo + 1) / 2;
        cpp_int p = pow_cpp_u32(mid, k);
        if (p <= val) lo = mid;
        else r = mid - 1;
    }
    return lo;
}

static inline uint32_t ceil_kth_root_cpp(const cpp_int& val, int k, uint32_t hi) {
    uint32_t x = floor_kth_root_cpp(val, k, hi);
    if (pow_cpp_u32(x, k) < val) x++;
    return x;
}

static py::list search_dfs_2sum_impl(
    int k,
    int m_terms,
    uint32_t min_b,
    uint32_t max_b,
    uint32_t max_val,
    uint32_t mod_M,
    uint32_t p1,
    uint32_t p2,
    uint32_t p3,
    uint32_t max_solutions,
    bool primitive_grouping,
    py::object progress_cb
) {
    if (k <= 0) throw std::invalid_argument("k must be positive");
    if (m_terms < 4) throw std::invalid_argument("m_terms must be >= 4 (needs at least a 2-sum completion)");
    if (min_b < 2) min_b = 2;
    if (max_b < min_b) return py::list();

    if (max_val == 0) max_val = max_b;
    if (max_val < 2) max_val = 2;

    if (mod_M != 65536u) throw std::invalid_argument("mod_M must be 65536 (2^16) in this backend");
    if (p1 >= 65536u || p2 >= 65536u || p3 >= 65536u) throw std::invalid_argument("p1,p2,p3 must be < 65536");
    if (max_val > 65535u) throw std::invalid_argument("DFS 2-sum backend supports max_val <= 65535");

    const uint64_t total_b = (uint64_t)max_b - (uint64_t)min_b + 1;

    // Precompute residues for 0..max_val
    std::vector<uint16_t> powM(max_val + 1);
    std::vector<uint16_t> powP1(max_val + 1);
    std::vector<uint16_t> powP2(max_val + 1);
    std::vector<uint16_t> powP3(max_val + 1);
    for (uint32_t i = 0; i <= max_val; ++i) {
        powM[i]  = (uint16_t)powmod_u32(i, k, mod_M);
        powP1[i] = (uint16_t)powmod_u32(i, k, p1);
        powP2[i] = (uint16_t)powmod_u32(i, k, p2);
        powP3[i] = (uint16_t)powmod_u32(i, k, p3);
    }

    const bool use_u128 = can_use_u128(max_val, max_b, k, m_terms);

    std::vector<u128> powExact128;
    std::vector<cpp_int> powExactCpp;
    if (use_u128) {
        powExact128.resize(max_val + 1);
        for (uint32_t i = 0; i <= max_val; ++i) powExact128[i] = pow_u128_fast(i, k);
    } else {
        powExactCpp.resize(max_val + 1);
        for (uint32_t i = 0; i <= max_val; ++i) powExactCpp[i] = pow_cpp_u32(i, k);
    }

    // Build 2-sum table: all pairs (a<=b) up to max_val, keyed by residues.
    // Memory O(n^2) and sorted for fast lookup + max-term cutoff.
    const uint64_t pair_count_est = ((uint64_t)max_val * (uint64_t)(max_val + 1)) / 2;
    std::vector<PairEntry> pairs;
    pairs.reserve((size_t)std::min<uint64_t>(pair_count_est, 5'000'000ull)); // reserve hint

    for (uint32_t a = 1; a <= max_val; ++a) {
        const uint16_t aM = powM[a], a1 = powP1[a], a2 = powP2[a], a3 = powP3[a];
        for (uint32_t b = a; b <= max_val; ++b) {
            const uint16_t sM = (uint16_t)(aM + powM[b]);

            uint32_t s1 = (uint32_t)a1 + (uint32_t)powP1[b]; if (s1 >= p1) s1 -= p1;
            uint32_t s2 = (uint32_t)a2 + (uint32_t)powP2[b]; if (s2 >= p2) s2 -= p2;
            uint32_t s3 = (uint32_t)a3 + (uint32_t)powP3[b]; if (s3 >= p3) s3 -= p3;

            const uint64_t key = pack_key_4(sM, (uint16_t)s1, (uint16_t)s2, (uint16_t)s3);
            pairs.push_back(PairEntry{key, (uint16_t)a, (uint16_t)b, (uint16_t)b});
        }
    }

    std::sort(pairs.begin(), pairs.end(), [](const PairEntry& x, const PairEntry& y) {
        if (x.key != y.key) return x.key < y.key;
        if (x.maxab != y.maxab) return x.maxab < y.maxab;
        if (x.b != y.b) return x.b < y.b;
        return x.a < y.a;
    });

    auto range_for_key = [&](uint64_t key) {
        auto lo = std::lower_bound(pairs.begin(), pairs.end(), key, [](const PairEntry& e, uint64_t k){ return e.key < k; });
        auto hi = std::upper_bound(pairs.begin(), pairs.end(), key, [](uint64_t k, const PairEntry& e){ return k < e.key; });
        return std::make_pair(lo, hi);
    };

    py::list out;
    std::unordered_set<uint64_t> seen;
    if (primitive_grouping) seen.reserve((size_t)max_solutions * 2);

    // Helper: emit solution
    auto emit_solution = [&](std::vector<uint32_t>& terms_desc, uint32_t B) {
        // Convert to ascending (canonical)
        std::vector<uint32_t> terms = terms_desc;
        std::sort(terms.begin(), terms.end());

        uint32_t outB = B;
        if (primitive_grouping) {
            uint32_t g = gcd_terms_and_b(terms, B);
            if (g > 1) {
                for (auto& v : terms) v /= g;
                outB /= g;
            }
            uint64_t sig = signature_terms_b(terms, outB);
            if (seen.find(sig) != seen.end()) return;
            seen.insert(sig);
        }

        py::tuple lhsT((size_t)m_terms);
        for (int i = 0; i < m_terms; ++i) lhsT[i] = terms[i];
        out.append(py::make_tuple(lhsT, outB));
    };

    // Iterate b
    uint64_t done = 0;
    for (uint32_t B = min_b; B <= max_b; ++B) {
        done++;
        progress_call(progress_cb, done, total_b);

        const uint32_t rhs_limit = (B < max_val ? B : max_val);

        const uint16_t tM  = (uint16_t)powmod_u32(B, k, mod_M);
        const uint16_t t1  = (uint16_t)powmod_u32(B, k, p1);
        const uint16_t t2  = (uint16_t)powmod_u32(B, k, p2);
        const uint16_t t3  = (uint16_t)powmod_u32(B, k, p3);

        // Target exact
        u128 target128 = 0;
        cpp_int targetCpp = 0;
        if (use_u128) target128 = pow_u128_fast(B, k);
        else targetCpp = pow_cpp_u32(B, k);

        // DFS over nonincreasing terms (descending), with 2-sum completion
        std::vector<uint32_t> current;
        current.reserve((size_t)m_terms);

        if (use_u128) {
            auto dfs = [&](auto&& self,
                           uint32_t prev,
                           u128 rem,
                           uint16_t sM, uint16_t s1, uint16_t s2, uint16_t s3,
                           int remaining) -> void {
                if ((uint32_t)out.size() >= max_solutions) return;

                if (remaining == 2) {
                    // Need residues for the final pair
                    const uint16_t needM = (uint16_t)(tM - sM);
                    const uint16_t need1 = sub_mod_u16(t1, s1, p1);
                    const uint16_t need2 = sub_mod_u16(t2, s2, p2);
                    const uint16_t need3 = sub_mod_u16(t3, s3, p3);
                    const uint64_t key = pack_key_4(needM, need1, need2, need3);

                    auto range = range_for_key(key);
                     auto lo = range.first;
                     auto hi = range.second;if (lo == hi) return;

                    const uint32_t limit = std::min(prev, rhs_limit);
                    for (auto it = lo; it != hi; ++it) {
                        if (it->maxab > limit) break; // sorted by maxab
                        const uint32_t a = it->a;
                        const uint32_t b = it->b;
                        u128 s = powExact128[a] + powExact128[b];
                        if (s != rem) continue;

                        // Build full term list (still descending if we append b then a)
                        current.push_back(b);
                        current.push_back(a);
                        emit_solution(current, B);
                        current.pop_back();
                        current.pop_back();

                        if ((uint32_t)out.size() >= max_solutions) return;
                    }
                    return;
                }

                // Minimal possible sum with remaining terms is remaining * 1^k = remaining
                if (rem < (u128)remaining) return;

                uint32_t ub0 = std::min(prev, rhs_limit);
                uint32_t ub = std::min(ub0, floor_kth_root_u128(rem, k, ub0));

                // Lower bound: need remaining*x^k >= rem  => x >= ceil_root(ceil(rem/remaining), k)
                u128 div = (rem + (u128)remaining - 1) / (u128)remaining;
                uint32_t lb = ceil_kth_root_u128(div, k, ub);
                if (lb < 1) lb = 1;
                if (lb > ub) return;

                for (uint32_t x = ub; ; --x) {
                    u128 px = powExact128[x];
                    if (px <= rem) {
                        current.push_back(x);

                        uint16_t nM = (uint16_t)(sM + powM[x]);
                        uint32_t nn1 = (uint32_t)s1 + (uint32_t)powP1[x]; if (nn1 >= p1) nn1 -= p1;
                        uint32_t nn2 = (uint32_t)s2 + (uint32_t)powP2[x]; if (nn2 >= p2) nn2 -= p2;
                        uint32_t nn3 = (uint32_t)s3 + (uint32_t)powP3[x]; if (nn3 >= p3) nn3 -= p3;

                        self(self, x, rem - px, nM, (uint16_t)nn1, (uint16_t)nn2, (uint16_t)nn3, remaining - 1);

                        current.pop_back();
                        if ((uint32_t)out.size() >= max_solutions) return;
                    }
                    if (x == lb) break;
                }
            };

            dfs(dfs, rhs_limit, target128, 0, 0, 0, 0, m_terms);
        } else {
            // cpp_int variant (slower, but supports very large B^k)
            auto dfs = [&](auto&& self,
                           uint32_t prev,
                           cpp_int& rem,
                           uint16_t sM, uint16_t s1, uint16_t s2, uint16_t s3,
                           int remaining) -> void {
                if ((uint32_t)out.size() >= max_solutions) return;

                if (remaining == 2) {
                    const uint16_t needM = (uint16_t)(tM - sM);
                    const uint16_t need1 = sub_mod_u16(t1, s1, p1);
                    const uint16_t need2 = sub_mod_u16(t2, s2, p2);
                    const uint16_t need3 = sub_mod_u16(t3, s3, p3);
                    const uint64_t key = pack_key_4(needM, need1, need2, need3);

                    auto range = range_for_key(key);
                     auto lo = range.first;
                     auto hi = range.second;if (lo == hi) return;

                    const uint32_t limit = std::min(prev, rhs_limit);
                    for (auto it = lo; it != hi; ++it) {
                        if (it->maxab > limit) break;
                        const uint32_t a = it->a;
                        const uint32_t b = it->b;
                        cpp_int s = powExactCpp[a] + powExactCpp[b];
                        if (s != rem) continue;

                        current.push_back(b);
                        current.push_back(a);
                        emit_solution(current, B);
                        current.pop_back();
                        current.pop_back();

                        if ((uint32_t)out.size() >= max_solutions) return;
                    }
                    return;
                }

                // rem must be >= remaining
                if (rem < remaining) return;

                uint32_t ub0 = std::min(prev, rhs_limit);
                uint32_t ub = std::min(ub0, floor_kth_root_cpp(rem, k, ub0));

                // div = ceil(rem/remaining)
                cpp_int div = (rem + (remaining - 1)) / remaining;
                uint32_t lb = ceil_kth_root_cpp(div, k, ub);
                if (lb < 1) lb = 1;
                if (lb > ub) return;

                for (uint32_t x = ub; ; --x) {
                    const cpp_int& px = powExactCpp[x];
                    if (px <= rem) {
                        current.push_back(x);

                        uint16_t nM = (uint16_t)(sM + powM[x]);
                        uint32_t nn1 = (uint32_t)s1 + (uint32_t)powP1[x]; if (nn1 >= p1) nn1 -= p1;
                        uint32_t nn2 = (uint32_t)s2 + (uint32_t)powP2[x]; if (nn2 >= p2) nn2 -= p2;
                        uint32_t nn3 = (uint32_t)s3 + (uint32_t)powP3[x]; if (nn3 >= p3) nn3 -= p3;

                        rem -= px;
                        self(self, x, rem, nM, (uint16_t)nn1, (uint16_t)nn2, (uint16_t)nn3, remaining - 1);
                        rem += px;

                        current.pop_back();
                        if ((uint32_t)out.size() >= max_solutions) return;
                    }
                    if (x == lb) break;
                }
            };

            cpp_int rem = targetCpp;
            dfs(dfs, rhs_limit, rem, 0, 0, 0, 0, m_terms);
        }

        if ((uint32_t)out.size() >= max_solutions) break;
    }

    return out;
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

    m.def(
        "search_dfs_2sum",
        &search_dfs_2sum_impl,
        py::arg("k"),
        py::arg("m_terms"),
        py::arg("min_b"),
        py::arg("max_b"),
        py::arg("max_val") = 0,
        py::arg("mod_M") = 65536,
        py::arg("p1") = 65521,
        py::arg("p2") = 65519,
        py::arg("p3") = 65497,
        py::arg("max_solutions") = 50,
        py::arg("primitive_grouping") = false,
        py::arg("progress_cb") = py::none()
    );
}