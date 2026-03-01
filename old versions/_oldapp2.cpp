#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <boost/multiprecision/cpp_int.hpp>
#include <boost/multiprecision/cpp_int.hpp>
#include <boost/multiprecision/cpp_int.hpp>

#include <algorithm>
#include <cstdint>
#include <numeric>
#include <stdexcept>
#include <vector>

namespace py = pybind11;
using boost::multiprecision::cpp_int;
using u128 = boost::multiprecision::uint128_t;

// -------------------- helpers --------------------

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

static inline void progress_call(const py::object& cb, uint64_t done, uint64_t total) {
    if (cb.is_none()) return;
    py::gil_scoped_acquire gil;
    cb((uint64_t)done, (uint64_t)total);
}

// Decide whether exact values fit safely in u128:
// Need: max_b^k fits, and m*max_val^k fits.
static inline bool can_use_u128(uint32_t max_val, uint32_t max_b, int k, int m) {
    cpp_int limit = cpp_int(1);
    limit <<= 128;
    cpp_int a = pow_cpp_u32(max_val, k);
    cpp_int b = pow_cpp_u32(max_b, k);
    cpp_int mm = m;
    return (a < limit) && (b < limit) && (a * mm < limit);
}

// -------------------- dynamic bitset for sieve primes --------------------

struct Bits {
    uint32_t nbits = 0;
    std::vector<uint64_t> w;

    Bits() = default;
    explicit Bits(uint32_t bits): nbits(bits), w((bits + 63u) / 64u, 0ull) {}

    inline void set(uint32_t i) {
        w[i >> 6] |= (1ull << (i & 63));
    }
    inline bool test(uint32_t i) const {
        return (w[i >> 6] >> (i & 63)) & 1ull;
    }
};

// rotate-left by shift within [0,nbits)
static inline Bits rotl_bits(const Bits& b, uint32_t shift) {
    if (b.nbits == 0) return Bits();
    shift %= b.nbits;
    if (shift == 0) return b;

    Bits out(b.nbits);
    const uint32_t words = (uint32_t)b.w.size();
    const uint32_t word_shift = shift >> 6;
    const uint32_t bit_shift  = shift & 63;

    for (uint32_t i = 0; i < words; ++i) {
        uint64_t lo = b.w[i];
        uint32_t j1 = (i + word_shift) % words;
        uint32_t j2 = (j1 + 1) % words;

        if (bit_shift == 0) {
            out.w[j1] |= lo;
        } else {
            out.w[j1] |= (lo << bit_shift);
            out.w[j2] |= (lo >> (64 - bit_shift));
        }
    }

    // mask off extra bits in last word
    uint32_t extra = (words * 64u) - b.nbits;
    if (extra) {
        out.w[words - 1] &= (~0ull) >> extra;
    }
    return out;
}

static inline void or_inplace(Bits& dst, const Bits& src) {
    if (dst.w.size() != src.w.size()) return;
    for (size_t i = 0; i < dst.w.size(); ++i) dst.w[i] |= src.w[i];
}

// For a given prime p, precompute reachable[t] bitsets for t=0..m where
// reachable[t][r] == 1 iff r is representable as sum of t kth-power residues mod p.
struct PrimeSieveDP {
    uint32_t p = 0;
    std::vector<uint32_t> residues;        // distinct x^k mod p
    std::vector<Bits> reachable;           // size m+1, each Bits(p)
    std::vector<uint32_t> term_residue;    // term_residue[x] = x^k mod p for x in 0..max_val
};

// -------------------- 2-sum bucket table --------------------

struct PairEntry {
    uint16_t a;
    uint16_t b;
    uint16_t r1; // sum mod p1
    uint16_t r2; // sum mod p2
};

// Build buckets indexed by (sum mod 65536). Each entry stores (a,b) and residues mod p1/p2.
// Terms are stored as a<=b.
static inline void build_pair_buckets(
    uint32_t limit,
    const std::vector<uint16_t>& powM,
    const std::vector<uint16_t>& powP1,
    const std::vector<uint16_t>& powP2,
    uint32_t p1,
    uint32_t p2,
    std::vector<std::vector<PairEntry>>& buckets
) {
    buckets.assign(65536u, {});
    // rough reserve: total pairs ~ limit*(limit+1)/2 spread across 65536 buckets, so small per bucket
    for (uint32_t b = 1; b <= limit; ++b) {
        uint16_t bM  = powM[b];
        uint16_t bP1 = powP1[b];
        uint16_t bP2 = powP2[b];
        for (uint32_t a = 1; a <= b; ++a) {
            uint16_t sumM = (uint16_t)(powM[a] + bM); // wrap mod 2^16
            uint32_t tmp1 = (uint32_t)powP1[a] + (uint32_t)bP1;
            uint16_t sum1 = (uint16_t)(tmp1 % p1);
            uint32_t tmp2 = (uint32_t)powP2[a] + (uint32_t)bP2;
            uint16_t sum2 = (uint16_t)(tmp2 % p2);
            buckets[sumM].push_back(PairEntry{(uint16_t)a,(uint16_t)b,sum1,sum2});
        }
    }
}

// -------------------- main solver --------------------

static py::list search_dfs_2sum_impl(
    int k,
    int m,
    uint32_t min_b,
    uint32_t max_b,
    uint32_t max_val,
    uint32_t mod_M,
    uint32_t p1,
    uint32_t p2,
    const std::vector<uint32_t>& sieve_primes,
    uint32_t max_solutions,
    bool primitive_grouping,
    py::object progress_cb
) {
    if (k <= 0) throw std::invalid_argument("k must be positive");
    if (m < 3) throw std::invalid_argument("m must be >= 3");
    if (min_b < 2) min_b = 2;
    if (max_b < min_b) return py::list();

    if (max_val == 0) max_val = max_b; // match Python semantics (terms <= b unless capped)
    if (max_val < 2) max_val = 2;
    if (mod_M != 65536u) throw std::invalid_argument("mod_M must be 65536 (2^16) for this backend");
    if (p1 == 0 || p2 == 0 || p1 >= 65536u || p2 >= 65536u) {
        throw std::invalid_argument("p1 and p2 must be in 1..65535");
    }

    const uint64_t total_b = (uint64_t)max_b - (uint64_t)min_b + 1;

    // Precompute residues for 0..max_val
    std::vector<uint16_t> powM(max_val + 1), powP1(max_val + 1), powP2(max_val + 1);
    for (uint32_t i = 0; i <= max_val; ++i) {
        powM[i]  = (uint16_t)powmod_u32(i, k, mod_M);
        powP1[i] = (uint16_t)powmod_u32(i, k, p1);
        powP2[i] = (uint16_t)powmod_u32(i, k, p2);
    }

    // Precompute sieve DPs (optional; helps prune DFS)
    std::vector<PrimeSieveDP> sieveDPs;
    sieveDPs.reserve(sieve_primes.size());
    for (uint32_t p : sieve_primes) {
        if (p < 2) continue;
        // Skip huge p to avoid heavy DP; correctness unaffected (just less pruning)
        if (p > 4096u) continue;

        PrimeSieveDP dp;
        dp.p = p;
        dp.term_residue.assign(max_val + 1, 0u);

        std::vector<uint8_t> seen(p, 0);
        for (uint32_t x = 0; x < p; ++x) {
            uint32_t r = powmod_u32(x, k, p);
            if (!seen[r]) { seen[r] = 1; dp.residues.push_back(r); }
        }
        std::sort(dp.residues.begin(), dp.residues.end());

        for (uint32_t x = 0; x <= max_val; ++x) dp.term_residue[x] = powmod_u32(x, k, p);

        dp.reachable.resize((size_t)m + 1);
        dp.reachable[0] = Bits(p);
        dp.reachable[0].set(0);

        for (int t = 1; t <= m; ++t) {
            Bits acc(p);
            const Bits& prev = dp.reachable[t - 1];
            for (uint32_t r : dp.residues) {
                Bits shifted = rotl_bits(prev, r);
                or_inplace(acc, shifted);
            }
            dp.reachable[t] = std::move(acc);
        }

        sieveDPs.push_back(std::move(dp));
    }

    const bool use_u128 = can_use_u128(max_val, max_b, k, m);

    std::vector<u128> powU;
    std::vector<u128> powBU;
    std::vector<cpp_int> powBig;
    std::vector<cpp_int> powBBig;

    if (use_u128) {
        powU.resize(max_val + 1);
        for (uint32_t i = 0; i <= max_val; ++i) powU[i] = pow_u128_fast(i, k);
        powBU.resize(max_b + 1);
        for (uint32_t i = 0; i <= max_b; ++i) powBU[i] = pow_u128_fast(i, k);
    } else {
        powBig.resize(max_val + 1);
        for (uint32_t i = 0; i <= max_val; ++i) powBig[i] = pow_cpp_u32(i, k);
        powBBig.resize(max_b + 1);
        for (uint32_t i = 0; i <= max_b; ++i) powBBig[i] = pow_cpp_u32(i, k);
    }

    // Pair buckets: build once at max_val.
    // NOTE: For correctness when b < max_val, we will filter pairs by b<=limit during lookup.
    std::vector<std::vector<PairEntry>> pairBuckets;
    build_pair_buckets(max_val, powM, powP1, powP2, p1, p2, pairBuckets);

    // Storage for results (dedup)
    py::list out;
    std::vector<uint64_t> sigs; sigs.reserve(max_solutions);

    auto have_sig = [&](uint64_t s) -> bool {
        return std::find(sigs.begin(), sigs.end(), s) != sigs.end();
    };

    auto push_solution = [&](std::vector<uint32_t> terms, uint32_t B) {
        std::sort(terms.begin(), terms.end()); // ensure nondecreasing

        if (primitive_grouping) {
            uint32_t g = B;
            for (uint32_t x : terms) g = gcd_u32(g, x);
            if (g > 1) {
                for (auto& x : terms) x /= g;
                B /= g;
            }
        }

        // signature for dedup
        uint64_t h = 1469598103934665603ull;
        for (uint32_t v : terms) {
            h ^= (uint64_t)v + 0x9e3779b97f4a7c15ull + (h<<6) + (h>>2);
        }
        h ^= (uint64_t)B + 0x9e3779b97f4a7c15ull + (h<<6) + (h>>2);

        if (have_sig(h)) return;
        sigs.push_back(h);

        py::list lhs;
        for (uint32_t v : terms) lhs.append((int)v);
        out.append(py::make_tuple(lhs, (int)B));
    };

    // Main loop over B
    for (uint32_t B = min_b; B <= max_b; ++B) {
        const uint32_t rhs_limit = std::min<uint32_t>(max_val, B);

        // Target residues
        const uint16_t targetM  = (uint16_t)powmod_u32(B, k, mod_M);
        const uint16_t targetP1 = (uint16_t)powmod_u32(B, k, p1);
        const uint16_t targetP2 = (uint16_t)powmod_u32(B, k, p2);

        // Target exact
        u128 targetU = 0;
        cpp_int targetBig = 0;
        if (use_u128) targetU = powBU[B];
        else targetBig = powBBig[B];

        // DFS arrays
        std::vector<uint32_t> chosen;
        chosen.reserve(m);

        // For sieve primes, keep partial residues per DP
        std::vector<uint32_t> partialSieve;
        partialSieve.assign(sieveDPs.size(), 0u);

        // Lambda to test sieve-prime feasibility for remaining terms count r_total
        auto sieve_ok = [&](int r_total, const std::vector<uint32_t>& partialRes) -> bool {
            for (size_t i = 0; i < sieveDPs.size(); ++i) {
                const auto& dp = sieveDPs[i];
                uint32_t need = (uint32_t)powmod_u32(B, k, dp.p);
                uint32_t pr = partialRes[i] % dp.p;
                need = (need + dp.p - pr) % dp.p;
                if (r_total >= 0 && r_total <= m) {
                    if (!dp.reachable[(size_t)r_total].test(need)) return false;
                }
            }
            return true;
        };

        // DFS recursion (choose m-2 terms, finish with 2-sum)
        std::function<void(int depth_left, uint32_t max_term,
                           u128 sumU, uint16_t sumM, uint16_t sumP1, uint16_t sumP2,
                           cpp_int* sumBig)> dfs;

        dfs = [&](int depth_left, uint32_t max_term,
                  u128 sumU, uint16_t sumM, uint16_t sumP1, uint16_t sumP2,
                  cpp_int* sumBigPtr) {
            if ((uint32_t)out.size() >= max_solutions) return;
            if (max_term < 1) return;

            const int remaining_total = depth_left + 2; // includes final pair
            if (!sieveDPs.empty()) {
                if (!sieve_ok(remaining_total, partialSieve)) return;
            }

            // remaining exact
            if (use_u128) {
                if (sumU > targetU) return;
                u128 rem = targetU - sumU;

                // upper bound: remaining_total * max_term^k must be >= rem
                u128 maxPow = powU[std::min<uint32_t>(max_term, rhs_limit)];
                u128 ub = maxPow * (u128)remaining_total;
                if (ub < rem) return;

                if (depth_left == 0) {
                    // Need a pair (a,b) with a<=b<=max_term and exact sum=rem and matching residues
                    uint16_t needM  = (uint16_t)(targetM - sumM); // wraps mod 2^16
                    uint16_t need1  = (uint16_t)((uint32_t)targetP1 + p1 - sumP1); if (need1 >= p1) need1 = (uint16_t)(need1 % p1);
                    uint16_t need2  = (uint16_t)((uint32_t)targetP2 + p2 - sumP2); if (need2 >= p2) need2 = (uint16_t)(need2 % p2);

                    auto& bucket = pairBuckets[needM];
                    for (const auto& pe : bucket) {
                        if (pe.b > max_term) continue;      // enforce nondecreasing
                        if (pe.b > rhs_limit) continue;     // enforce term<=b
                        if (pe.r1 != need1) continue;
                        if (pe.r2 != need2) continue;

                        u128 s2 = powU[pe.a] + powU[pe.b];
                        if (s2 == rem) {
                            std::vector<uint32_t> terms = chosen;
                            terms.push_back(pe.a);
                            terms.push_back(pe.b);

                            // verify exact (paranoia)
                            u128 check = 0;
                            for (uint32_t x : terms) check += powU[x];
                            if (check == targetU) {
                                push_solution(std::move(terms), B);
                                if ((uint32_t)out.size() >= max_solutions) return;
                            }
                        }
                    }
                    return;
                }

                // Choose next term t (descending) up to max_term and rhs_limit, also t^k <= rem.
                uint32_t hi = std::min<uint32_t>(max_term, rhs_limit);
                while (hi > 1 && powU[hi] > rem) --hi;

                for (uint32_t t = hi; t >= 1; --t) {
                    u128 newSum = sumU + powU[t];
                    if (newSum > targetU) continue;

                    chosen.push_back(t);
                    // update sieve partial residues
                    for (size_t i = 0; i < sieveDPs.size(); ++i) {
                        const auto& dp = sieveDPs[i];
                        partialSieve[i] += dp.term_residue[t];
                        partialSieve[i] %= dp.p;
                    }

                    dfs(depth_left - 1, t, newSum,
                        (uint16_t)(sumM + powM[t]),
                        (uint16_t)((sumP1 + powP1[t]) % p1),
                        (uint16_t)((sumP2 + powP2[t]) % p2),
                        nullptr);

                    // revert
                    for (size_t i = 0; i < sieveDPs.size(); ++i) {
                        const auto& dp = sieveDPs[i];
                        // undo (mod p)
                        uint32_t v = dp.term_residue[t] % dp.p;
                        partialSieve[i] = (partialSieve[i] + dp.p - v) % dp.p;
                    }
                    chosen.pop_back();

                    if (t == 1) break;
                }
            } else {
                // big-int path
                cpp_int& sumBig = *sumBigPtr;
                if (sumBig > targetBig) return;
                cpp_int rem = targetBig - sumBig;

                cpp_int maxPow = powBig[std::min<uint32_t>(max_term, rhs_limit)];
                cpp_int ub = maxPow * remaining_total;
                if (ub < rem) return;

                if (depth_left == 0) {
                    uint16_t needM  = (uint16_t)(targetM - sumM);
                    uint16_t need1  = (uint16_t)((uint32_t)targetP1 + p1 - sumP1); if (need1 >= p1) need1 = (uint16_t)(need1 % p1);
                    uint16_t need2  = (uint16_t)((uint32_t)targetP2 + p2 - sumP2); if (need2 >= p2) need2 = (uint16_t)(need2 % p2);

                    auto& bucket = pairBuckets[needM];
                    for (const auto& pe : bucket) {
                        if (pe.b > max_term) continue;
                        if (pe.b > rhs_limit) continue;
                        if (pe.r1 != need1) continue;
                        if (pe.r2 != need2) continue;

                        cpp_int s2 = powBig[pe.a] + powBig[pe.b];
                        if (s2 == rem) {
                            std::vector<uint32_t> terms = chosen;
                            terms.push_back(pe.a);
                            terms.push_back(pe.b);

                            // verify exact
                            cpp_int check = 0;
                            for (uint32_t x : terms) check += powBig[x];
                            if (check == targetBig) {
                                push_solution(std::move(terms), B);
                                if ((uint32_t)out.size() >= max_solutions) return;
                            }
                        }
                    }
                    return;
                }

                uint32_t hi = std::min<uint32_t>(max_term, rhs_limit);
                while (hi > 1 && powBig[hi] > rem) --hi;

                for (uint32_t t = hi; t >= 1; --t) {
                    cpp_int newSum = sumBig + powBig[t];
                    if (newSum > targetBig) continue;

                    chosen.push_back(t);
                    for (size_t i = 0; i < sieveDPs.size(); ++i) {
                        const auto& dp = sieveDPs[i];
                        partialSieve[i] += dp.term_residue[t];
                        partialSieve[i] %= dp.p;
                    }

                    cpp_int saved = sumBig;
                    sumBig = newSum;
                    dfs(depth_left - 1, t, 0, // sumU unused
                        (uint16_t)(sumM + powM[t]),
                        (uint16_t)((sumP1 + powP1[t]) % p1),
                        (uint16_t)((sumP2 + powP2[t]) % p2),
                        &sumBig);
                    sumBig = saved;

                    for (size_t i = 0; i < sieveDPs.size(); ++i) {
                        const auto& dp = sieveDPs[i];
                        uint32_t v = dp.term_residue[t] % dp.p;
                        partialSieve[i] = (partialSieve[i] + dp.p - v) % dp.p;
                    }
                    chosen.pop_back();

                    if (t == 1) break;
                }
            }
        };

        if (use_u128) {
            dfs(m - 2, rhs_limit, 0, 0, 0, 0, nullptr);
        } else {
            cpp_int sum0 = 0;
            dfs(m - 2, rhs_limit, 0, 0, 0, 0, &sum0);
        }

        progress_call(progress_cb, (uint64_t)(B - min_b + 1), total_b);
        if ((uint32_t)out.size() >= max_solutions) break;
    }

    return out;
}

PYBIND11_MODULE(euler_ext, m) {
    m.doc() = "Euler sum-of-powers search (RAM-safe DFS + 2-sum completion + sieve pruning)";

    m.def(
        "search_dfs_2sum",
        &search_dfs_2sum_impl,
        py::arg("k"),
        py::arg("m"),
        py::arg("min_b"),
        py::arg("max_b"),
        py::arg("max_val") = 0,
        py::arg("mod_M") = 65536,
        py::arg("p1") = 65521,
        py::arg("p2") = 65519,
        py::arg("sieve_primes") = std::vector<uint32_t>{},
        py::arg("max_solutions") = 50,
        py::arg("primitive_grouping") = false,
        py::arg("progress_cb") = py::none()
    );
}
