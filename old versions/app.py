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

def sieve_b_values(k: int, m: int, primes: List[int], N: int) -> List[int]:
    """Keep b where b^k mod p is representable as sum of m k-th power residues mod p for all primes p."""
    tables = []
    for p in primes:
        res = kth_residues_mask(k, p)
        reachable = mfold_sumset_mask(res, m, p)
        tables.append((p, reachable))

    good = []
    for b in range(2, N + 1):
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
    """Number of nondecreasing r-tuples in [1..N] equals C(N+r-1, r)."""
    return math.comb(N + r - 1, r)

def precompute_pows(k: int, N: int) -> List[int]:
    """pows[i] = i^k (Python int is arbitrary precision)."""
    return [0] + [pow(i, k) for i in range(1, N + 1)]

def enumerate_sums(
    pows: List[int],
    length: int,
    max_val: int,
    start: int = 1,
    prefix: Optional[List[int]] = None,
    running_sum: int = 0,
    target: Optional[int] = None,
) -> Iterable[Tuple[Tuple[int, ...], int]]:
    """
    Enumerate nondecreasing tuples of given length from [start..max_val],
    yielding (tuple, sum(pows[x])).

    target is an upper bound for the sum (used for pruning).
    Fast-path: if prefix is empty and length in {1,2,3,4}, uses iterative loops.
    """
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
    """
    Build a table mapping sum -> witness of length r.

    Optimizations:
      - packed witnesses (r<=4, max_val<=65535)
      - optional bucketing by sum % bucket_mod to reduce dict size and cache misses
    """
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
# Solvers: brute force (validation) + MITM (research)
# ============================================================

def brute_solutions(k: int, m: int, N: int, max_solutions: int = 50) -> List[Tuple[Tuple[int, ...], int]]:
    """
    Brute-force find solutions to a1^k+...+am^k=b^k with nondecreasing a's and a_i <= b <= N.
    Only feasible for small N/m.
    """
    pows = precompute_pows(k, N)
    sols: List[Tuple[Tuple[int, ...], int]] = []
    t0 = time.perf_counter()

    for b in range(2, N + 1):
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
    N: int,
    *,
    split_r: Optional[int] = None,
    primes: Optional[List[int]] = None,
    bucket_mod: Optional[int] = None,
    store_all: bool = False,
    max_solutions: int = 50,
    progress: bool = True,
    primitive_grouping: bool = False,
) -> List[Tuple[Tuple[int, ...], int]]:
    """
    Meet-in-the-middle for a1^k+...+am^k=b^k with nondecreasing tuples and a_i <= b <= N.

    Optimizations included:
      - forced good splits for m=4 (r=2) and m=8 (r=4)
      - optional sieve on b
      - optional bucketing (bucket_mod)
      - packed witnesses in sum table (r<=4, N<=65535)
      - iterative enumeration for lengths 1..4
      - dedupe of exact duplicates, and optional primitive grouping
    """
    if m < 2:
        raise ValueError("m must be >= 2")

    # Upgrade A: force good splits for common cases
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

    pows = precompute_pows(k, N)

    # Optional sieve on b
    if primes:
        b_candidates = sieve_b_values(k, m, primes, N)
        if progress:
            print(f"[mitm] sieve kept {len(b_candidates):,}/{max(0, N-1):,} b-candidates using primes={primes}")
    else:
        b_candidates = list(range(2, N + 1))

    # Build r-sum table once
    table_info = build_sum_table(pows, r, N, bucket_mod=bucket_mod, store_all=store_all, progress=progress)

    sols: List[Tuple[Tuple[int, ...], int]] = []
    seen = set()
    t0 = time.perf_counter()

    for idx, b in enumerate(b_candidates, 1):
        target = pows[b]

        # Enumerate (m-r)-sums up to b (enforces a_i <= b)
        for rhs_tup, rhs_sum in enumerate_sums(pows, m - r, b, target=target):
            need = target - rhs_sum
            if need < 0:
                continue

            hit = table_lookup(table_info, need)
            if hit is None:
                continue

            def add_solution(lhs_tuple: Tuple[int, ...]) -> bool:
                full = tuple(sorted(lhs_tuple + rhs_tup))
                key = normalize_solution(full, b) if primitive_grouping else (full, b)
                if key in seen:
                    return False
                seen.add(key)
                sols.append((full, b))
                return True

            if store_all and isinstance(hit, list):
                for w in hit:
                    lhs_tuple = materialize_witness(table_info, w)
                    add_solution(lhs_tuple)
                    if len(sols) >= max_solutions:
                        if progress:
                            dt = time.perf_counter() - t0
                            print(f"[mitm] reached max_solutions={max_solutions} time={dt:.2f}s")
                        return sols
            else:
                lhs_tuple = materialize_witness(table_info, hit)  # type: ignore
                add_solution(lhs_tuple)
                if len(sols) >= max_solutions:
                    if progress:
                        dt = time.perf_counter() - t0
                        print(f"[mitm] reached max_solutions={max_solutions} time={dt:.2f}s")
                    return sols

        if progress and (idx % max(1, len(b_candidates) // 10) == 0):
            print(f"[mitm] progress: {idx:,}/{len(b_candidates):,} b checked, sols={len(sols)}")

    if progress:
        dt = time.perf_counter() - t0
        print(f"[mitm] solutions={len(sols)} time={dt:.2f}s")

    return sols

# ============================================================
# Console UI: Validation + Research modes
# ============================================================

def normalize_solution_list(sols: List[Tuple[Tuple[int, ...], int]]) -> set:
    return set(sols)

def run_validation_mode():
    print("\n=== VALIDATION MODE ===")
    print("Goal: brute-force vs MITM on SMALL bounds and assert they match.\n")

    k = int(input("Exponent k (e.g., 2, 4, 5): ").strip())
    m = int(input("LHS term count m (Euler-style use m=k-1): ").strip())
    N = int(input("Max value N (small! e.g. 20..60): ").strip())

    brute_est = estimate_nondec_tuple_count(N, m) * max(1, N - 1)
    print(f"\n[info] rough brute tuple count scale ~ {brute_est:,} (very rough)")
    if brute_est > 20_000_000:
        print("[warning] This brute run may be slow. Consider smaller N or smaller m.\n")

    split = input(f"MITM split r (blank for default floor(m/2)={m//2}): ").strip()
    split_r = int(split) if split else None

    max_s = input("Max solutions to collect (default 50): ").strip()
    max_solutions = int(max_s) if max_s else 50

    print("\nRunning brute force...")
    brute = brute_solutions(k, m, N, max_solutions=max_solutions)

    print("\nRunning meet-in-the-middle...")
    mitm = mitm_solutions(
        k, m, N,
        split_r=split_r,
        primes=None,
        bucket_mod=None,
        store_all=True,        # collect all witnesses in validation for exact comparison
        max_solutions=max_solutions,
        progress=True,
        primitive_grouping=False,
    )

    brute_set = normalize_solution_list(brute)
    mitm_set = normalize_solution_list(mitm)

    print("\n=== VALIDATION RESULT ===")
    if brute_set == mitm_set:
        print("PASS! brute and MITM matched exactly on this parameter set.")
    else:
        only_brute = brute_set - mitm_set
        only_mitm = mitm_set - brute_set
        print("FAIL! mismatch detected.")
        if only_brute:
            print(f"  Found by brute only ({len(only_brute)}): {list(only_brute)[:5]}")
        if only_mitm:
            print(f"  Found by MITM only ({len(only_mitm)}): {list(only_mitm)[:5]}")

    if brute:
        print("\nExample solutions (lhs -> b):")
        for lhs, b in brute[:min(10, len(brute))]:
            print(f"  {lhs} -> {b}")

def run_research_mode():
    print("\n=== RESEARCH MODE ===")
    print("Goal: MITM search with optional b-sieving, bucketing, packed tables, and primitive/multiple flags.\n")

    k = int(input("Exponent k (e.g., 5 or 9): ").strip())
    m = int(input("LHS term count m (Euler-style uses m=k-1): ").strip())
    N = int(input("Max value N (start small for big k; e.g., 60): ").strip())

    split = input(f"MITM split r (blank for auto; default floor(m/2)={m//2}): ").strip()
    split_r = int(split) if split else None

    primes_str = input("Sieve primes (comma list) or blank for none (e.g., 19,37,73,109): ").strip()
    primes = [int(x) for x in primes_str.split(",") if x.strip()] if primes_str else None

    bucket_str = input("Bucket modulus (blank for auto; examples: 4096, 8192, 1155): ").strip()
    bucket_mod = int(bucket_str) if bucket_str else None

    # Defaults: force good splits, and enable bucket_mod for hard cases
    auto_r = (split_r if split_r is not None else (2 if m == 4 else 4 if m == 8 else (m // 2)))
    if split_r is None:
        split_r = auto_r

    if bucket_mod is None:
        # Upgrade B: default bucketing for hard regimes
        if k >= 7 or split_r >= 3 or N >= 200:
            bucket_mod = 4096

    store_all = input("Store ALL witnesses per sum? (y/N): ").strip().lower() == "y"
    grouping = input("Group scaled multiples into primitive classes? (y/N): ").strip().lower() == "y"

    max_s = input("Max solutions to collect (default 20): ").strip()
    max_solutions = int(max_s) if max_s else 20

    r = split_r
    est_table = estimate_nondec_tuple_count(N, r)
    est_other = estimate_nondec_tuple_count(N, m - r)
    print(f"\n[info] using split r={r} and (m-r)={m-r}")
    print(f"[info] estimated r-sum tuples ≈ {est_table:,}")
    print(f"[info] estimated (m-r)-sum tuples per b up to N ≈ {est_other:,}")
    print(f"[info] bucket_mod={'none' if bucket_mod is None else bucket_mod}")
    print("[note] Total work is roughly r-sums + (b loop)*(m-r sums). Start small for k>=7.\n")

    sols = mitm_solutions(
        k, m, N,
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
        print("No solutions found within the chosen bound.")
        return

    print(f"Found {len(sols)} solutions (showing up to 20):")
    for lhs, b in sols[:20]:
        g = solution_gcd(lhs, b)
        tag = "[primitive]" if g == 1 else f"[multiple ×{g}]"
        print(f"  {lhs} -> {b}   {tag}")

def main():
    while True:
        print("\n==============================")
        print(" Euler / Sum-of-Powers Search ")
        print("==============================")
        print("1) Validation mode (brute vs MITM)")
        print("2) Research mode (MITM + sieve/buckets + packed tables)")
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
            print("Closing Program...")
            return
        else:
            print("Invalid choice.")

if __name__ == "__main__":
    main()