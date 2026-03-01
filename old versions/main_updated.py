# main_updated.py
# Speed-optimized (without sacrificing mathematical legitimacy) sum-of-powers search.
#
# Key upgrades vs your current main.py:
#   1) Range control on b: min_b..max_b, AND we only precompute up to max_b (real speed win).
#   2) Hard-regime fast path: for m=8, split r=4 -> precompute RHS 4-sums once and reuse for all b.
#   3) Better progress stats: per-b scan rate + periodic global ETA.
#   4) Keeps all mathematical correctness: we do NOT restrict a_i ranges (only b range if you choose).
#
# Compatible with Python 3.8+.

import time
import math
from typing import Dict, List, Tuple, Iterable, Optional, Union

# ============================================================
# Solution classification (primitive vs scaled)
# ============================================================

def solution_gcd(lhs: Tuple[int, ...], b: int) -> int:
    g = b
    for x in lhs:
        g = math.gcd(g, x)
    return g

def normalize_solution(lhs: Tuple[int, ...], b: int) -> Tuple[Tuple[int, ...], int]:
    g = solution_gcd(lhs, b)
    if g == 1:
        return lhs, b
    return tuple(x // g for x in lhs), b // g

# ============================================================
# Modular sieve helpers (fast b-candidate filtering)
# ============================================================

def kth_residues_mask(k: int, p: int) -> int:
    mask = 0
    for x in range(p):
        mask |= 1 << pow(x, k, p)
    return mask

def add_cyclic(mask: int, shift: int, p: int) -> int:
    allbits = (1 << p) - 1
    shift %= p
    return ((mask << shift) | (mask >> (p - shift))) & allbits

def mfold_sumset_mask(res_mask: int, m: int, p: int) -> int:
    reachable = 1  # {0}
    for _ in range(m):
        new = 0
        rm = res_mask
        while rm:
            lsb = rm & -rm
            r = (lsb.bit_length() - 1)
            new |= add_cyclic(reachable, r, p)
            rm -= lsb
        reachable = new
    return reachable

def sieve_b_values(k: int, m: int, primes: List[int], max_b: int) -> List[int]:
    """Keep b where b^k mod p is representable as sum of m k-th power residues mod p for all primes p."""
    tables = []
    for p in primes:
        res = kth_residues_mask(k, p)
        reachable = mfold_sumset_mask(res, m, p)
        tables.append((p, reachable))

    good: List[int] = []
    for b in range(2, max_b + 1):
        ok = True
        for p, reachable in tables:
            t = pow(b, k, p)
            if ((reachable >> t) & 1) == 0:
                ok = False
                break
        if ok:
            good.append(b)
    return good

# ============================================================
# Combinatorics + optimizations
#   - iterative enumeration for lengths 1..4 (much faster than recursion)
#   - packed witnesses for r<=4 (cuts memory/hash overhead)
# ============================================================

def estimate_nondec_tuple_count(N: int, r: int) -> int:
    return math.comb(N + r - 1, r)

def precompute_pows(k: int, max_b: int) -> List[int]:
    return [0] + [pow(i, k) for i in range(1, max_b + 1)]

def enumerate_sums(
    pows: List[int],
    length: int,
    max_val: int,
    start: int = 1,
    prefix: Optional[List[int]] = None,
    running_sum: int = 0,
    target: Optional[int] = None,
) -> Iterable[Tuple[Tuple[int, ...], int]]:
    if prefix is None:
        prefix = []

    if target is not None and running_sum > target:
        return

    if length == 0:
        yield (tuple(prefix), running_sum)
        return

    # lower bound prune
    if target is not None and start <= max_val:
        min_possible = running_sum + length * pows[start]
        if min_possible > target:
            return

    # fast-path (only when no prefix has been built yet)
    if len(prefix) == 0 and length in (1, 2, 3, 4):
        if length == 1:
            for a in range(start, max_val + 1):
                s = running_sum + pows[a]
                if target is not None and s > target:
                    break
                yield ((a,), s)
            return

        if length == 2:
            for a in range(start, max_val + 1):
                sa = running_sum + pows[a]
                if target is not None and sa + pows[a] > target:
                    break
                for b in range(a, max_val + 1):
                    s = sa + pows[b]
                    if target is not None and s > target:
                        break
                    yield ((a, b), s)
            return

        if length == 3:
            for a in range(start, max_val + 1):
                sa = running_sum + pows[a]
                if target is not None and sa + 2 * pows[a] > target:
                    break
                for b in range(a, max_val + 1):
                    sab = sa + pows[b]
                    if target is not None and sab + pows[b] > target:
                        break
                    for c in range(b, max_val + 1):
                        s = sab + pows[c]
                        if target is not None and s > target:
                            break
                        yield ((a, b, c), s)
            return

        # length == 4
        for a in range(start, max_val + 1):
            sa = running_sum + pows[a]
            if target is not None and sa + 3 * pows[a] > target:
                break
            for b in range(a, max_val + 1):
                sab = sa + pows[b]
                if target is not None and sab + 2 * pows[b] > target:
                    break
                for c in range(b, max_val + 1):
                    sabc = sab + pows[c]
                    if target is not None and sabc + pows[c] > target:
                        break
                    for d in range(c, max_val + 1):
                        s = sabc + pows[d]
                        if target is not None and s > target:
                            break
                        yield ((a, b, c, d), s)
        return

    # general recursive fallback
    for x in range(start, max_val + 1):
        s2 = running_sum + pows[x]
        if target is not None and s2 > target:
            break
        prefix.append(x)
        for out in enumerate_sums(
            pows=pows,
            length=length - 1,
            max_val=max_val,
            start=x,
            prefix=prefix,
            running_sum=s2,
            target=target,
        ):
            yield out
        prefix.pop()

# ---- Packed witness (r <= 4, max_val <= 65535) ----

_PACK_SHIFT = 16
_PACK_MASK = (1 << _PACK_SHIFT) - 1

def can_pack_witness(r: int, max_val: int) -> bool:
    return (1 <= r <= 4) and (max_val <= _PACK_MASK)

def pack_witness(tup: Tuple[int, ...]) -> int:
    x = 0
    for v in tup:
        x = (x << _PACK_SHIFT) | (v & _PACK_MASK)
    return x

def unpack_witness(x: int, r: int) -> Tuple[int, ...]:
    out = [0] * r
    for i in range(r - 1, -1, -1):
        out[i] = x & _PACK_MASK
        x >>= _PACK_SHIFT
    return tuple(out)

# ============================================================
# Sum table builder (optional bucketing + packed witnesses)
# ============================================================

Witness = Tuple[int, ...]
PackedWitness = int
TableValue = Union[Witness, PackedWitness, List[Witness], List[PackedWitness]]
SumTableBuckets = List[Dict[int, TableValue]]
SumTablePlain = Dict[int, TableValue]
SumTable = Union[SumTablePlain, SumTableBuckets]

class SumTableInfo:
    def __init__(self, table: SumTable, r: int, bucket_mod: Optional[int], packed: bool):
        self.table = table
        self.r = r
        self.bucket_mod = bucket_mod
        self.packed = packed

def build_sum_table(
    pows: List[int],
    r: int,
    max_val: int,
    *,
    bucket_mod: Optional[int] = None,
    store_all: bool = False,
    progress: bool = True,
) -> SumTableInfo:
    if r <= 0:
        raise ValueError("r must be >= 1 for build_sum_table")

    packed = can_pack_witness(r, max_val)
    t0 = time.perf_counter()
    est = estimate_nondec_tuple_count(max_val, r)

    if progress:
        bm = "none" if bucket_mod is None else str(bucket_mod)
        store = "packed" if packed else "tuples"
        print(f"[build_sum_table] r={r}, max_val={max_val}, est tuples≈{est:,}, bucket_mod={bm}, store={store}")

    def to_store(t: Tuple[int, ...]) -> Union[Tuple[int, ...], int]:
        return pack_witness(t) if packed else t

    if bucket_mod is None:
        table: Dict[int, TableValue] = {}
        for tup, s in enumerate_sums(pows, r, max_val):
            w = to_store(tup)
            if store_all:
                if s not in table:
                    table[s] = [w]  # type: ignore
                else:
                    table[s].append(w)  # type: ignore
            else:
                table.setdefault(s, w)
        if progress:
            dt = time.perf_counter() - t0
            print(f"[build_sum_table] unique sums={len(table):,} time={dt:.2f}s")
        return SumTableInfo(table, r=r, bucket_mod=None, packed=packed)

    if bucket_mod <= 0:
        raise ValueError("bucket_mod must be positive")

    buckets: List[Dict[int, TableValue]] = [dict() for _ in range(bucket_mod)]
    for tup, s in enumerate_sums(pows, r, max_val):
        w = to_store(tup)
        b = s % bucket_mod
        d = buckets[b]
        if store_all:
            if s not in d:
                d[s] = [w]  # type: ignore
            else:
                d[s].append(w)  # type: ignore
        else:
            d.setdefault(s, w)

    if progress:
        dt = time.perf_counter() - t0
        nonempty = sum(1 for d in buckets if d)
        uniques = sum(len(d) for d in buckets)
        print(f"[build_sum_table] buckets={bucket_mod} nonempty={nonempty:,} unique sums={uniques:,} time={dt:.2f}s")

    return SumTableInfo(buckets, r=r, bucket_mod=bucket_mod, packed=packed)

def table_lookup(info: SumTableInfo, need: int) -> Optional[TableValue]:
    if info.bucket_mod is None:
        return info.table.get(need)  # type: ignore
    b = need % info.bucket_mod
    return info.table[b].get(need)  # type: ignore

def materialize_witness(info: SumTableInfo, w: Union[Witness, PackedWitness]) -> Witness:
    if info.packed:
        return unpack_witness(int(w), info.r)
    return w  # type: ignore

# ============================================================
# Fast path for m=8,r=4: precompute RHS 4-sums once
# ============================================================

RHSItem = Tuple[int, int, Tuple[int, int, int, int]]  # (max_elem, sum, tuple)

def precompute_rhs_4sums(pows: List[int], max_val: int, progress: bool = True) -> Tuple[List[int], List[RHSItem]]:
    """
    Precompute all nondecreasing 4-tuples (a<=b<=c<=d<=max_val) and their sums.
    Returns lists sorted by max_elem so for each b we can scan only prefix where max_elem<=b.
    """
    t0 = time.perf_counter()
    items: List[RHSItem] = []
    for tup, s in enumerate_sums(pows, 4, max_val):
        items.append((tup[3], s, (tup[0], tup[1], tup[2], tup[3])))
    items.sort(key=lambda x: x[0])
    max_elems = [it[0] for it in items]
    if progress:
        dt = time.perf_counter() - t0
        print(f"[rhs] precomputed 4-sums: count={len(items):,} time={dt:.2f}s")
    return max_elems, items

def bisect_right(a: List[int], x: int) -> int:
    lo, hi = 0, len(a)
    while lo < hi:
        mid = (lo + hi) // 2
        if a[mid] <= x:
            lo = mid + 1
        else:
            hi = mid
    return lo

# ============================================================
# Solvers
# ============================================================

def brute_solutions(k: int, m: int, max_b: int, max_solutions: int = 50) -> List[Tuple[Tuple[int, ...], int]]:
    pows = precompute_pows(k, max_b)
    sols: List[Tuple[Tuple[int, ...], int]] = []
    t0 = time.perf_counter()

    for b in range(2, max_b + 1):
        target = pows[b]
        for lhs, s in enumerate_sums(pows, m, b, target=target):
            if s == target:
                sols.append((lhs, b))
                if len(sols) >= max_solutions:
                    dt = time.perf_counter() - t0
                    print(f"[brute] hit max_solutions={max_solutions} after {dt:.2f}s")
                    return sols

    dt = time.perf_counter() - t0
    print(f"[brute] solutions={len(sols)} time={dt:.2f}s")
    return sols

def mitm_solutions(
    k: int,
    m: int,
    max_b: int,
    *,
    min_b: int = 2,
    split_r: Optional[int] = None,
    primes: Optional[List[int]] = None,
    bucket_mod: Optional[int] = None,
    store_all: bool = False,
    max_solutions: int = 50,
    progress: bool = True,
    primitive_grouping: bool = False,
) -> List[Tuple[Tuple[int, ...], int]]:
    """
    Meet-in-the-middle for a1^k+...+am^k=b^k with nondecreasing tuples and a_i <= b.

    Mathematical legitimacy:
      - We only restrict b to [min_b..max_b] (user-controlled).
      - We do NOT restrict a_i beyond the required a_i <= b and nondecreasing ordering.

    Speed upgrades:
      - Precompute pows up to max_b
      - Build LHS table only up to max_b (not any larger N)
      - For m=8,r=4: precompute RHS 4-sums once and reuse for all b
      - Better progress stats (scan rates + ETA)
    """
    if m < 2:
        raise ValueError("m must be >= 2")
    if min_b < 2:
        min_b = 2
    if max_b < min_b:
        raise ValueError("max_b must be >= min_b")

    # Force good splits
    if split_r is None:
        if m == 4:
            split_r = 2
        elif m == 8:
            split_r = 4
        else:
            split_r = m // 2

    r = split_r
    if not (1 <= r <= m - 1):
        raise ValueError("split_r must satisfy 1 <= r <= m-1")

    # Precompute up to max_b (this is the big legitimacy-preserving speed win)
    pows = precompute_pows(k, max_b)

    # Candidate b list (with sieve, then apply min_b)
    if primes:
        b_candidates = sieve_b_values(k, m, primes, max_b)
        if progress:
            print(f"[mitm] sieve kept {len(b_candidates):,}/{max(0, max_b-1):,} b-candidates using primes={primes}")
    else:
        b_candidates = list(range(2, max_b + 1))

    if min_b > 2:
        before = len(b_candidates)
        b_candidates = [b for b in b_candidates if b >= min_b]
        if progress:
            print(f"[mitm] applied min_b={min_b}: kept {len(b_candidates):,}/{before:,} candidates")

    # Build LHS r-sum table once, only up to max_b
    table_info = build_sum_table(pows, r, max_b, bucket_mod=bucket_mod, store_all=store_all, progress=progress)

    # Optional RHS precompute for m=8,r=4
    rhs_max_elems: Optional[List[int]] = None
    rhs_items: Optional[List[RHSItem]] = None
    if m == 8 and r == 4:
        rhs_max_elems, rhs_items = precompute_rhs_4sums(pows, max_b, progress=progress)

    sols: List[Tuple[Tuple[int, ...], int]] = []
    seen = set()

    # Stats
    t_all0 = time.perf_counter()
    rhs_scanned_total = 0
    b_done = 0
    last_global_print = t_all0

    for idx, b in enumerate(b_candidates, 1):
        b_done = idx
        t_b0 = time.perf_counter()
        target = pows[b]
        rhs_scanned = 0
        sols_before = len(sols)

        def add_solution(lhs_tuple: Tuple[int, ...], rhs_tuple: Tuple[int, ...]) -> None:
            full = tuple(sorted(lhs_tuple + rhs_tuple))
            key = normalize_solution(full, b) if primitive_grouping else (full, b)
            if key in seen:
                return
            seen.add(key)
            sols.append((full, b))

        if rhs_items is not None and rhs_max_elems is not None:
            limit = bisect_right(rhs_max_elems, b)
            for j in range(limit):
                _, rhs_sum, rhs_tup4 = rhs_items[j]
                rhs_scanned += 1
                need = target - rhs_sum
                if need < 0:
                    continue
                hit = table_lookup(table_info, need)
                if hit is None:
                    continue

                if store_all and isinstance(hit, list):
                    for w in hit:
                        lhs = materialize_witness(table_info, w)
                        add_solution(lhs, rhs_tup4)
                        if len(sols) >= max_solutions:
                            break
                else:
                    lhs = materialize_witness(table_info, hit)  # type: ignore
                    add_solution(lhs, rhs_tup4)

                if len(sols) >= max_solutions:
                    break
        else:
            # General path: enumerate RHS per b
            for rhs_tup, rhs_sum in enumerate_sums(pows, m - r, b, target=target):
                rhs_scanned += 1
                need = target - rhs_sum
                if need < 0:
                    continue
                hit = table_lookup(table_info, need)
                if hit is None:
                    continue

                if store_all and isinstance(hit, list):
                    for w in hit:
                        lhs = materialize_witness(table_info, w)
                        add_solution(lhs, rhs_tup)
                        if len(sols) >= max_solutions:
                            break
                else:
                    lhs = materialize_witness(table_info, hit)  # type: ignore
                    add_solution(lhs, rhs_tup)

                if len(sols) >= max_solutions:
                    break

        rhs_scanned_total += rhs_scanned

        # Per-b progress line (compact but informative)
        if progress:
            dt_b = time.perf_counter() - t_b0
            new_sols = len(sols) - sols_before
            rate = (rhs_scanned / dt_b) if dt_b > 1e-9 else 0.0
            print(f"[b={b}] rhs_scanned={rhs_scanned:,} rate={rate:,.0f}/s new_sols={new_sols} total_sols={len(sols)}")

        if len(sols) >= max_solutions:
            break

        # Global progress line every ~2 seconds or every 10% of candidates
        now = time.perf_counter()
        if progress and (now - last_global_print >= 2.0 or idx % max(1, len(b_candidates)//10) == 0):
            last_global_print = now
            dt_all = now - t_all0
            b_rate = (idx / dt_all) if dt_all > 1e-9 else 0.0
            rhs_rate = (rhs_scanned_total / dt_all) if dt_all > 1e-9 else 0.0

            remaining = len(b_candidates) - idx
            eta = (remaining / b_rate) if b_rate > 1e-9 else float("inf")
            eta_str = f"{eta:.0f}s" if eta != float("inf") else "∞"

            print(f"[mitm] {idx:,}/{len(b_candidates):,} b | b_rate={b_rate:.2f}/s | rhs_rate={rhs_rate:,.0f}/s | ETA~{eta_str} | sols={len(sols)}")

    if progress:
        dt_all = time.perf_counter() - t_all0
        b_rate = (b_done / dt_all) if dt_all > 1e-9 else 0.0
        rhs_rate = (rhs_scanned_total / dt_all) if dt_all > 1e-9 else 0.0
        print(f"[mitm] done: b_checked={b_done:,} time={dt_all:.2f}s b_rate={b_rate:.2f}/s rhs_scanned={rhs_scanned_total:,} rhs_rate={rhs_rate:,.0f}/s sols={len(sols)}")

    return sols

# ============================================================
# Console UI
# ============================================================

def normalize_solution_list(sols: List[Tuple[Tuple[int, ...], int]]) -> set:
    return set(sols)

def run_validation_mode() -> None:
    print("\n=== VALIDATION MODE ===")
    print("Goal: brute-force vs MITM on SMALL bounds and assert they match.\n")

    k = int(input("Exponent k (e.g., 2, 4, 5): ").strip())
    m = int(input("LHS term count m (Euler-style use m=k-1): ").strip())
    max_b = int(input("Max value max_b (small! e.g. 20..60): ").strip())

    brute_est = estimate_nondec_tuple_count(max_b, m) * max(1, max_b - 1)
    print(f"\n[info] rough brute tuple count scale ~ {brute_est:,} (very rough)")
    if brute_est > 20_000_000:
        print("[warning] This brute run may be slow. Consider smaller max_b or smaller m.\n")

    split = input(f"MITM split r (blank for auto; default floor(m/2)={m//2}): ").strip()
    split_r = int(split) if split else None

    max_s = input("Max solutions to collect (default 50): ").strip()
    max_solutions = int(max_s) if max_s else 50

    print("\nRunning brute force...")
    brute = brute_solutions(k, m, max_b, max_solutions=max_solutions)

    print("\nRunning meet-in-the-middle...")
    mitm = mitm_solutions(
        k, m, max_b,
        min_b=2,
        split_r=split_r,
        primes=None,
        bucket_mod=None,
        store_all=True,
        max_solutions=max_solutions,
        progress=True,
        primitive_grouping=False,
    )

    brute_set = normalize_solution_list(brute)
    mitm_set = normalize_solution_list(mitm)

    print("\n=== VALIDATION RESULT ===")
    if brute_set == mitm_set:
        print("PASS ✅  brute and MITM matched exactly on this parameter set.")
    else:
        only_brute = brute_set - mitm_set
        only_mitm = mitm_set - brute_set
        print("FAIL ❌  mismatch detected.")
        if only_brute:
            print(f"  Found by brute only ({len(only_brute)}): {list(only_brute)[:5]}")
        if only_mitm:
            print(f"  Found by MITM only ({len(only_mitm)}): {list(only_mitm)[:5]}")

def run_research_mode() -> None:
    print("\n=== RESEARCH MODE ===")
    print("Goal: fast MITM search with optional b-sieving, bucketing, b-range control, and strong progress stats.\n")

    k = int(input("Exponent k (e.g., 5 or 9): ").strip())
    m = int(input("LHS term count m (Euler-style uses m=k-1): ").strip())

    min_b_str = input("Minimum b to start from (blank for 2): ").strip()
    min_b = int(min_b_str) if min_b_str else 2

    max_b = int(input("Maximum b (max_b): ").strip())
    if max_b < min_b:
        print("max_b must be >= min_b.")
        return

    split = input(f"MITM split r (blank for auto; default floor(m/2)={m//2}, forces 2 for m=4 and 4 for m=8): ").strip()
    split_r = int(split) if split else None

    primes_str = input("Sieve primes (comma list) or blank for none (recommended for k>=7). Example: 19,37,73,109: ").strip()
    primes = [int(x) for x in primes_str.split(",") if x.strip()] if primes_str else None

    bucket_str = input("Bucket modulus (blank for auto; examples: 4096, 8192, 16384): ").strip()
    bucket_mod = int(bucket_str) if bucket_str else None

    # Auto defaults:
    if split_r is None:
        split_r = 4 if m == 8 else (2 if m == 4 else (m // 2))

    if bucket_mod is None:
        # Cache-friendly powers of two. Larger bucket is often better for 4-sums.
        if (m == 8 and split_r == 4 and (k >= 7 or max_b >= 200)):
            bucket_mod = 16384
        elif (k >= 7 or split_r >= 3 or max_b >= 200):
            bucket_mod = 4096
        else:
            bucket_mod = None

    store_all = input("Store ALL witnesses per sum? (y/N): ").strip().lower() == "y"
    grouping = input("Group scaled multiples into primitive classes? (y/N): ").strip().lower() == "y"

    max_s = input("Max solutions to collect (default 20): ").strip()
    max_solutions = int(max_s) if max_s else 20

    r = split_r
    print("\n=== PARAMETERS ===")
    print(f"k={k}, m={m}, min_b={min_b}, max_b={max_b}, split_r={r}")
    print(f"bucket_mod={'none' if bucket_mod is None else bucket_mod}")
    print(f"sieve primes={'none' if not primes else primes}")
    print(f"store_all={store_all}, primitive_grouping={grouping}, max_solutions={max_solutions}")
    print("\n[note] Mathematical legitimacy is preserved: we only restrict b-range, not a_i.\n")

    sols = mitm_solutions(
        k, m, max_b,
        min_b=min_b,
        split_r=r,
        primes=primes,
        bucket_mod=bucket_mod,
        store_all=store_all,
        max_solutions=max_solutions,
        progress=True,
        primitive_grouping=grouping,
    )

    print("\n=== RESEARCH OUTPUT ===")
    if not sols:
        print("No solutions found within the chosen b-range.")
        return

    print(f"Found {len(sols)} solutions (showing up to 20):")
    for lhs, b in sols[:20]:
        g = solution_gcd(lhs, b)
        tag = "[primitive]" if g == 1 else f"[multiple ×{g}]"
        print(f"  {lhs} -> {b}   {tag}")

def main() -> None:
    while True:
        print("\n==============================")
        print(" Euler / Sum-of-Powers Search ")
        print("==============================")
        print("1) Validation mode (brute vs MITM)")
        print("2) Research mode (fast MITM + b-range + strong progress stats)")
        print("3) Quick known check: 27^5+84^5+110^5+133^5 == 144^5")
        print("4) Exit")
        choice = input("Select: ").strip()

        if choice == "1":
            run_validation_mode()
        elif choice == "2":
            run_research_mode()
        elif choice == "3":
            lhs = pow(27, 5) + pow(84, 5) + pow(110, 5) + pow(133, 5)
            rhs = pow(144, 5)
            print("Check:", lhs == rhs)
            print("LHS:", lhs)
            print("RHS:", rhs)
        elif choice == "4":
            print("Bye.")
            return
        else:
            print("Invalid choice.")

if __name__ == "__main__":
    main()
