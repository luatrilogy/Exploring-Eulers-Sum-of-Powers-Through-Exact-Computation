
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <boost/multiprecision/cpp_int.hpp>

#include <cstdint>
#include <vector>
#include <algorithm>
#include <stdexcept>
#include <unordered_set>
#include <string>

namespace py = pybind11;
using boost::multiprecision::cpp_int;

// -------------------- utilities --------------------

static inline uint32_t powmod_u32(uint32_t a, int k, uint32_t mod) {
    uint64_t res = 1;
    uint64_t base = a % mod;
    int e = k;
    while (e > 0) {
        if (e & 1) res = (res * base) % mod;
        base = (base * base) % mod;
        e >>= 1;
    }
    return (uint32_t)res;
}

static inline cpp_int pow_cpp(uint32_t a, int k) {
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

struct Witness4 {
    uint32_t a, b, c, d;
};

struct Entry {
    // fingerprint packs (modP1, modP2, modM) into 48 bits:
    // fp = (p1<<32) | (p2<<16) | m
    uint64_t fp;
    Witness4 w;
};

static inline uint64_t pack_fp(uint32_t sP1, uint32_t sP2, uint32_t sM) {
    return ((uint64_t)(sP1 & 0xFFFFu) << 32) | ((uint64_t)(sP2 & 0xFFFFu) << 16) | (uint64_t)(sM & 0xFFFFu);
}

static inline uint32_t fp_p1(uint64_t fp){ return (uint32_t)((fp >> 32) & 0xFFFFu); }
static inline uint32_t fp_p2(uint64_t fp){ return (uint32_t)((fp >> 16) & 0xFFFFu); }
static inline uint32_t fp_m (uint64_t fp){ return (uint32_t)( fp        & 0xFFFFu); }

static inline uint64_t pack8_primitive(const std::vector<uint32_t>& xs) {
    // Hash-like packing for primitive grouping: compute gcd then normalize.
    // We use a simple rolling hash into 64-bit; only used for deduping output, not for correctness.
    uint64_t h = 1469598103934665603ull;
    for (auto v: xs) {
        h ^= (uint64_t)v + 0x9e3779b97f4a7c15ull + (h<<6) + (h>>2);
    }
    return h;
}

static inline uint32_t gcd_u32(uint32_t a, uint32_t b) {
    while (b) { uint32_t t = a % b; a = b; b = t; }
    return a;
}

static inline uint32_t gcd_many(const std::vector<uint32_t>& xs) {
    uint32_t g = 0;
    for (auto v: xs) g = (g==0? v : gcd_u32(g, v));
    return g == 0 ? 1 : g;
}

// -------------------- core: sieved + blocked m=8 --------------------

static py::list search_m8_blocked_sieved(
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
    if (mod_M != 65536u) throw std::invalid_argument("Currently mod_M must be 65536 (2^16) for packed fingerprints");
    if (p1 >= 65536u || p2 >= 65536u) throw std::invalid_argument("p1/p2 must be < 65536");
    if (min_b < 2) min_b = 2;
    if (max_b < min_b) return py::list();
    if (max_val < 2) max_val = max_b;

    // Precompute power residues for 1..max_val
    std::vector<uint16_t> powM(max_val + 1);
    std::vector<uint16_t> powP1(max_val + 1);
    std::vector<uint16_t> powP2(max_val + 1);
    for (uint32_t i=0; i<=max_val; ++i) {
        powM[i]  = (uint16_t)powmod_u32(i, k, mod_M);
        powP1[i] = (uint16_t)powmod_u32(i, k, p1);
        powP2[i] = (uint16_t)powmod_u32(i, k, p2);
    }

    // Build LHS fingerprint table ONCE (depends only on max_val, k, moduli).
    // This is the dominant memory structure; we keep it compact and capped.
    std::vector<std::vector<Entry>> buckets(mod_M);
    buckets.shrink_to_fit();

    uint64_t total_entries = 0;
    for (uint32_t a=1; a<=max_val; ++a) {
        uint16_t aM = powM[a], a1 = powP1[a], a2 = powP2[a];
        for (uint32_t b=a; b<=max_val; ++b) {
            uint16_t abM = (uint16_t)(aM + powM[b]);
            uint16_t ab1 = (uint16_t)((a1 + powP1[b]) % p1);
            uint16_t ab2 = (uint16_t)((a2 + powP2[b]) % p2);
            for (uint32_t c=b; c<=max_val; ++c) {
                uint16_t abcM = (uint16_t)(abM + powM[c]);
                uint16_t abc1 = (uint16_t)((ab1 + powP1[c]) % p1);
                uint16_t abc2 = (uint16_t)((ab2 + powP2[c]) % p2);
                for (uint32_t d=c; d<=max_val; ++d) {
                    uint16_t sM  = (uint16_t)(abcM + powM[d]);     // mod 2^16 wrap ok
                    uint16_t sP1 = (uint16_t)((abc1 + powP1[d]) % p1);
                    uint16_t sP2 = (uint16_t)((abc2 + powP2[d]) % p2);

                    Entry e;
                    e.fp = pack_fp(sP1, sP2, sM);
                    e.w  = Witness4{a,b,c,d};
                    buckets[sM].push_back(e);
                    total_entries++;

                    if (total_entries >= max_entries_cap) {
                        throw std::runtime_error("LHS table exceeded max_entries_cap; reduce max_val or increase cap.");
                    }
                }
            }
        }
    }

    // Sort each bucket by fp for fast scanning
    for (auto &vec : buckets) {
        std::sort(vec.begin(), vec.end(), [](const Entry& x, const Entry& y){ return x.fp < y.fp; });
    }

    // Output solutions
    py::list out;
    std::unordered_set<uint64_t> seen; // for primitive grouping
    seen.reserve((size_t)max_solutions * 2);

    auto verify_exact = [&](const Witness4& L, const Witness4& R, uint32_t b)->bool{
        cpp_int lhs = pow_cpp(L.a,k) + pow_cpp(L.b,k) + pow_cpp(L.c,k) + pow_cpp(L.d,k)
                    + pow_cpp(R.a,k) + pow_cpp(R.b,k) + pow_cpp(R.c,k) + pow_cpp(R.d,k);
        cpp_int rhs = pow_cpp(b,k);
        return lhs == rhs;
    };

    // Scan b in blocks (currently mainly for progress + future extensibility)
    for (uint32_t bStart = min_b; bStart <= max_b; ) {
        uint32_t bEnd = bStart + block_b - 1;
        if (bEnd > max_b) bEnd = max_b;

        for (uint32_t b = bStart; b <= bEnd; ++b) {
            if (!progress_cb.is_none()) progress_cb(b);

            uint16_t tM  = (uint16_t)powmod_u32(b, k, mod_M);
            uint16_t tP1 = (uint16_t)powmod_u32(b, k, p1);
            uint16_t tP2 = (uint16_t)powmod_u32(b, k, p2);

            // Enumerate RHS 4-sums (compact residue arithmetic first)
            for (uint32_t a=1; a<=max_val; ++a) {
                uint16_t aM = powM[a], a1 = powP1[a], a2 = powP2[a];
                for (uint32_t bb=a; bb<=max_val; ++bb) {
                    uint16_t abM = (uint16_t)(aM + powM[bb]);
                    uint16_t ab1 = (uint16_t)((a1 + powP1[bb]) % p1);
                    uint16_t ab2 = (uint16_t)((a2 + powP2[bb]) % p2);
                    for (uint32_t c=bb; c<=max_val; ++c) {
                        uint16_t abcM = (uint16_t)(abM + powM[c]);
                        uint16_t abc1 = (uint16_t)((ab1 + powP1[c]) % p1);
                        uint16_t abc2 = (uint16_t)((ab2 + powP2[c]) % p2);
                        for (uint32_t d=c; d<=max_val; ++d) {
                            uint16_t rM  = (uint16_t)(abcM + powM[d]);
                            uint16_t rP1 = (uint16_t)((abc1 + powP1[d]) % p1);
                            uint16_t rP2 = (uint16_t)((abc2 + powP2[d]) % p2);

                            // needed = target - rhs (mod each modulus)
                            uint16_t needM  = (uint16_t)(tM - rM);
                            uint16_t needP1 = (uint16_t)((tP1 + p1 - rP1) % p1);
                            uint16_t needP2 = (uint16_t)((tP2 + p2 - rP2) % p2);

                            uint64_t need_fp_prefix = ((uint64_t)(needP1 & 0xFFFFu) << 32) | ((uint64_t)(needP2 & 0xFFFFu) << 16);
                            // scan bucket needM for matching p1/p2 (and needM implied by bucket)
                            const auto &vec = buckets[needM];
                            if (vec.empty()) continue;

                            // fp is ordered; we can lower_bound by full fp
                            uint64_t key_fp = need_fp_prefix | (uint64_t)(needM & 0xFFFFu);
                            auto it = std::lower_bound(vec.begin(), vec.end(), key_fp,
                                [](const Entry& e, uint64_t key){ return e.fp < key; });

                            for (; it != vec.end(); ++it) {
                                if (fp_m(it->fp) != needM) break;
                                if (fp_p1(it->fp) != needP1) break; // because p1 is highest field, order groups by p1
                                if (fp_p2(it->fp) != needP2) continue;

                                Witness4 L = it->w;
                                Witness4 R = Witness4{a, bb, c, d};

                                if (!verify_exact(L, R, b)) continue;

                                // Construct output tuple of 8 terms (sorted nondecreasing by construction within each 4-tuple; cross-order not guaranteed)
                                std::vector<uint32_t> terms = {L.a,L.b,L.c,L.d,R.a,R.b,R.c,R.d};
                                std::sort(terms.begin(), terms.end());

                                if (primitive_grouping) {
                                    uint32_t g = gcd_many(terms);
                                    if (g > 1) {
                                        for (auto &v: terms) v /= g;
                                    }
                                    uint64_t sig = pack8_primitive(terms) ^ (uint64_t)b;
                                    if (seen.find(sig) != seen.end()) continue;
                                    seen.insert(sig);
                                }

                                py::tuple lhs8(8);
                                for (int i=0;i<8;i++) lhs8[i] = terms[i];
                                out.append(py::make_tuple(lhs8, b));

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

PYBIND11_MODULE(euler_ext, m) {
    m.doc() = "Euler sum-of-powers MITM core (C++/pybind11)";

    m.def(
        "search_m8_blocked_sieved",
        &search_m8_blocked_sieved,
        py::arg("k"),
        py::arg("min_b"),
        py::arg("max_b"),
        py::arg("max_val") = 0,
        py::arg("mod_M") = 65536,
        py::arg("p1") = 65521,   // largest 16-bit prime
        py::arg("p2") = 65519,
        py::arg("block_b") = 50,
        py::arg("max_solutions") = 50,
        py::arg("max_entries_cap") = (uint64_t)25'000'000, // ~25M entries cap by default
        py::arg("primitive_grouping") = true,
        py::arg("progress_cb") = py::none()
    );
}
