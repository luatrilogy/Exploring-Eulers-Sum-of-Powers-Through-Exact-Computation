# path = .\.venv\Scripts\python.exe 

import sys
import time
import math
from tqdm import tqdm
from typing import Dict, List, Tuple, Iterable, Optional, Union

try:
    import euler_ext
    print("OK")
    HAS_CPP = True
except Exception as e:
    euler_ext = None 
    print("FAIL")
    HAS_CPP = False
    CPP_IMPORT_ERROR = repr(e)

from tqdm import tqdm
import sys

# --- progress bar globals (shared by C++ backends) ---
pbar = None
last_done = 0
total_b = 0

def on_progress(done, tot=None):
    """Progress callback used by the C++ extension.

    The C++ extension calls progress_cb(done, total).
    """
    global pbar, last_done, total_b
    done = int(done)
    if tot is not None:
        total_b = int(tot)

    if pbar is None:
        pbar = tqdm(
            total=total_b if total_b else None,
            desc="MITM search (b)",
            unit="b",
            dynamic_ncols=True,
            file=sys.stdout,
            leave=True,
        )

    if done > last_done:
        pbar.update(done - last_done)
        last_done = done

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
# ============================================================
# For m=8 split 4+4, LHS is 4 terms:

def residues_of_kth_powers(k: int, mod: int) -> List[int]:
    """All residues x^k mod mod for x in 0..mod-1 (deduped)."""
    seen = set()
    for x in range(mod):
        seen.add(pow(x, k, mod))
    return sorted(seen)

def rotate_bitset(bits: int, shift: int, mod: int) -> int:
    """Rotate-left a mod-bit bitset stored in a Python int."""
    shift %= mod
    if shift == 0:
        return bits
    mask = (1 << mod) - 1
    return ((bits << shift) | (bits >> (mod - shift))) & mask

def mfold_sumset_bitmask(residues: List[int], m: int, mod: int) -> int:
    """
    Bitmask of residues achievable as sum of m residues (mod mod).
    If bit r is 1, residue r is achievable.
    """
    reachable = 1  # residue 0
    for _ in range(m):
        new = 0
        for r in residues:
            new |= rotate_bitset(reachable, r, mod)
        reachable = new
    return reachable

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
    """Enumerate nondecreasing tuples of given length (fast for length<=4)."""
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

    for x in range(start, max_val + 1):
        s2 = running_sum + pows[x]
        if target is not None and s2 > target:
            break
        prefix.append(x)
        for out in enumerate_sums(pows, length - 1, max_val, x, prefix, s2, target):
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
    max_witnesses_per_sum: int = 0,
    progress: bool = True,
) -> SumTableInfo:
    packed = can_pack_witness(r, max_val)
    t0 = time.perf_counter()
    est = estimate_nondec_tuple_count(max_val, r)

    if progress:
        bm = "none" if bucket_mod is None else str(bucket_mod)
        store = "packed" if packed else "tuples"
        cap = "" if (not store_all or max_witnesses_per_sum == 0) else f", cap={max_witnesses_per_sum}/sum"
        print(f"[build_sum_table] r={r}, max_val={max_val}, est tuples≈{est:,}, bucket_mod={bm}, store={store}{cap}")

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
                    lst = table[s]  # type: ignore
                    if (max_witnesses_per_sum == 0) or (len(lst) < max_witnesses_per_sum):
                        lst.append(w)  # type: ignore
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
                lst = d[s]  # type: ignore
                if (max_witnesses_per_sum == 0) or (len(lst) < max_witnesses_per_sum):
                    lst.append(w)  # type: ignore
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
# RHS 4-sums grouped by max_elem and sorted by sum (LEGIT cutoff)
# ============================================================

# Each RHS entry is (rhs_sum, rhs_tuple4)
RHSPlainEntry = Tuple[int, Tuple[int, int, int, int]]
RHSPlainGrouped = List[List[RHSPlainEntry]]  # index by max_elem (0..max_val)

# Residue-aware RHS entries (rhs_sum, rhs_sum_mod_p, rhs_sum_mod_M, tuple4)
RHSResidueEntry = Tuple[int, int, int, Tuple[int, int, int, int]]
RHSResidueGrouped = List[List[RHSResidueEntry]]

def precompute_rhs_4sums_grouped(pows: List[int], max_val: int, progress: bool = True) -> RHSPlainGrouped:
    """
    For each d in 1..max_val, store all 4-tuples with max_elem==d, as (sum, tuple),
    and sort each list by sum ascending.

    Then for a given b with target=b^k, we scan:
      for d=1..b: for (sum, tup) in group[d]:
                    if sum>target: break   (LEGIT cutoff)
                    otherwise test hit.
    """
    t0 = time.perf_counter()
    groups: RHSPlainGrouped = [[] for _ in range(max_val + 1)]
    count = 0
    for tup, s in enumerate_sums(pows, 4, max_val):
        d = tup[3]
        groups[d].append((s, (tup[0], tup[1], tup[2], tup[3])))
        count += 1

    # sort each group by sum for early cutoff
    nonempty = 0
    for d in range(1, max_val + 1):
        if groups[d]:
            groups[d].sort(key=lambda x: x[0])
            nonempty += 1

    if progress:
        dt = time.perf_counter() - t0
        print(f"[rhs] grouped 4-sums: total={count:,} groups_nonempty={nonempty:,} time={dt:.2f}s")
    return groups

# Solver
def cpp_m8_solutions(
    k: int,
    max_b: int,
    *,
    min_b: int = 2,
    mod_M: int = 65536,
    max_val: int = 0,
    max_solutions: int = 50,
    primitive_grouping: bool = False,
    sieve_primes: Optional[List[int]] = None,
) -> list[tuple[tuple[int, ...], int]]:
    """C++ backend (RAM-safe): DFS + 2-sum completion for m=8."""
    if not HAS_CPP:
        raise RuntimeError(f"C++ backend not available: {CPP_IMPORT_ERROR}")

    global pbar, last_done, total_b
    pbar = None
    last_done = 0
    total_b = max_b - min_b + 1

    try:
        raw = euler_ext.search_dfs_2sum(
            k=int(k),
            m=8,
            min_b=int(min_b),
            max_b=int(max_b),
            max_val=int(max_val) if max_val else 0,
            mod_M=int(mod_M),
            p1=65521,
            p2=65519,
            sieve_primes=[int(x) for x in (sieve_primes or [])],
            max_solutions=int(max_solutions),
            primitive_grouping=bool(primitive_grouping),
            progress_cb=on_progress,
        )
    finally:
        if pbar is not None:
            pbar.close()
            pbar = None

    sols: list[tuple[tuple[int, ...], int]] = []
    seen = set()
    for lhs, b in raw:
        lhs = tuple(int(x) for x in lhs)
        b = int(b)
        key = normalize_solution(lhs, b) if primitive_grouping else (lhs, b)
        if key in seen:
            continue
        seen.add(key)
        sols.append((lhs, b))
    return sols

def cpp_m6_solutions(
    k: int,
    max_b: int,
    *,
    min_b: int = 2,
    max_val: int = 0,
    max_solutions: int = 50,
    primitive_grouping: bool = False,
    sieve_primes: Optional[List[int]] = None,
):
    """C++ backend (RAM-safe): DFS + 2-sum completion for m=6.

    Note: this is not as fast as a dedicated 3+3 MITM, but it will not allocate O(n^3) tables,
    and it stays exact.
    """
    if not HAS_CPP:
        raise RuntimeError(f"C++ Backend Not Available: {CPP_IMPORT_ERROR}")

    global pbar, last_done, total_b
    pbar = None
    last_done = 0
    total_b = max_b - min_b + 1

    try:
        raw = euler_ext.search_dfs_2sum(
            k=int(k),
            m=6,
            min_b=int(min_b),
            max_b=int(max_b),
            max_val=int(max_val) if max_val else 0,
            mod_M=65536,
            p1=65521,
            p2=65519,
            sieve_primes=[int(x) for x in (sieve_primes or [])],
            max_solutions=int(max_solutions),
            primitive_grouping=bool(primitive_grouping),
            progress_cb=on_progress,
        )
    finally:
        if pbar is not None:
            pbar.close()
            pbar = None

    sols = []
    seen = set()
    for lhs, b in raw:
        lhs = tuple(int(x) for x in lhs)
        b = int(b)
        key = normalize_solution(lhs, b) if primitive_grouping else (lhs, b)
        if key in seen:
            continue
        seen.add(key)
        sols.append((lhs, b))
    return sols

def cpp_m9_solutions(
    k: int,
    max_b: int,
    *,
    min_b: int = 2,
    max_val: int = 0,
    max_solutions: int = 50,
    primitive_grouping: bool = False,
    sieve_primes: Optional[List[int]] = None,
):
    """C++ backend (RAM-safe): DFS + 2-sum completion for m=9 (intended for k=9)."""
    if not HAS_CPP:
        raise RuntimeError(f"C++ Backend Not Available: {CPP_IMPORT_ERROR}")

    global pbar, last_done, total_b
    pbar = None
    last_done = 0
    total_b = max_b - min_b + 1

    try:
        raw = euler_ext.search_dfs_2sum(
            k=int(k),
            m=9,
            min_b=int(min_b),
            max_b=int(max_b),
            max_val=int(max_val) if max_val else 0,
            mod_M=65536,
            p1=65521,
            p2=65519,
            sieve_primes=[int(x) for x in (sieve_primes or [])],
            max_solutions=int(max_solutions),
            primitive_grouping=bool(primitive_grouping),
            progress_cb=on_progress,
        )
    finally:
        if pbar is not None:
            pbar.close()
            pbar = None

    sols = []
    seen = set()
    for lhs, b in raw:
        lhs = tuple(int(x) for x in lhs)
        b = int(b)
        key = normalize_solution(lhs, b) if primitive_grouping else (lhs, b)
        if key in seen:
            continue
        seen.add(key)
        sols.append((lhs, b))
    return sols




def precompute_rhs_4sums_grouped_with_residues(
    pows: List[int],
    max_val: int,
    mod_p: int,
    mod_M: int,
) -> RHSResidueGrouped:
    groups: RHSResidueGrouped = [[] for _ in range(max_val + 1)]
    for tup, s in enumerate_sums(pows, 4, max_val):
        d = tup[3]
        groups[d].append((s, s % mod_p, s % mod_M, (tup[0], tup[1], tup[2], tup[3])))

    for d in range(1, max_val + 1):
        if groups[d]:
            groups[d].sort(key=lambda x: x[0])
    return groups


def brute_solutions(k: int, m: int, N: int, max_solutions: int = 50) -> List[Tuple[Tuple[int, ...], int]]:
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
    per_b_print_every: int = 1,
) -> List[Tuple[Tuple[int, ...], int]]:
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

    pows = precompute_pows(k, max_b)

    # b candidates
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

    # LHS table
    table_info = build_sum_table(
        pows, r, max_b,
        bucket_mod=bucket_mod,
        store_all=store_all,
        max_witnesses_per_sum=(0 if store_all else 1),
        progress=progress,
    )

        # Residue filters (only used for the hard case m=8 split 4+4)
    mod_p: Optional[int] = None
    mod_M: Optional[int] = None
    lhs_res_p_mask: int = 0
    lhs_res_M_mask: int = 0

    if m == 8 and r == 4:
        mod_p = 73        # try 19 / 37 / 73
        mod_M = 4096      # try 2048 / 4096 / 8192 (power of two)

        lhs_res_p_mask = mfold_sumset_bitmask(
            residues_of_kth_powers(k, mod_p), r, mod_p
        )
        lhs_res_M_mask = mfold_sumset_bitmask(
            residues_of_kth_powers(k, mod_M), r, mod_M
        )

        if progress:
            print(f"[filter] LHS residues mod {mod_p}: {lhs_res_p_mask.bit_count()}/{mod_p}")
            print(f"[filter] LHS residues mod {mod_M}: {lhs_res_M_mask.bit_count()}/{mod_M}")

# RHS grouped cutoff for the key hard case m=8,r=4
    rhs_groups: Optional[RHSResidueGrouped] = None
    if m == 8 and r == 4:
        assert mod_p is not None and mod_M is not None
        rhs_groups = precompute_rhs_4sums_grouped_with_residues(pows, max_b, mod_p, mod_M)

    sols: List[Tuple[Tuple[int, ...], int]] = []
    seen = set()

    # stats
    t0 = time.perf_counter()
    rhs_scanned_total = 0
    last_global_print = t0
    b_done = 0

    for idx, b in enumerate(b_candidates, 1):
        b_done = idx
        target = pows[b]
        t_b0 = time.perf_counter()
        rhs_scanned = 0
        sols_before = len(sols)

        def add_solution(lhs_tuple: Tuple[int, ...], rhs_tuple: Tuple[int, ...]) -> None:
            full = tuple(sorted(lhs_tuple + rhs_tuple))
            key = normalize_solution(full, b) if primitive_grouping else (full, b)
            if key in seen:
                return
            seen.add(key)
            sols.append((full, b))

        if rhs_groups is not None:
            assert mod_p is not None and mod_M is not None
            t_res_p = target % mod_p
            t_res_M = target % mod_M

            for d in range(1, b + 1):
                grp = rhs_groups[d]
                if not grp:
                    continue

                for rhs_sum, rhs_rp, rhs_rM, rhs_tup4 in grp:
                    if rhs_sum > target:
                        break   

                    rhs_scanned += 1

                    # ---- residue filters ----
                    need_rp = (t_res_p - rhs_rp) % mod_p
                    if ((lhs_res_p_mask >> need_rp) & 1) == 0:
                        continue

                    need_rM = (t_res_M - rhs_rM) % mod_M
                    if ((lhs_res_M_mask >> need_rM) & 1) == 0:
                        continue
                    # -------------------------

                    need = target - rhs_sum
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

                if len(sols) >= max_solutions:
                    break
        else:
            # General path (no grouped cutoff)
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

        if progress and (per_b_print_every > 0) and (idx % per_b_print_every == 0):
            dt_b = time.perf_counter() - t_b0
            new_sols = len(sols) - sols_before
            rate = (rhs_scanned / dt_b) if dt_b > 1e-9 else 0.0
            print(f"[b={b}] rhs_scanned={rhs_scanned:,} rate={rate:,.0f}/s new_sols={new_sols} total_sols={len(sols)}")

        if len(sols) >= max_solutions:
            break

        now = time.perf_counter()
        if progress and (now - last_global_print >= 2.0 or idx % max(1, len(b_candidates)//per_b_print_every) == 0):
            last_global_print = now
            dt = now - t0
            b_rate = (idx / dt) if dt > 1e-9 else 0.0
            rhs_rate = (rhs_scanned_total / dt) if dt > 1e-9 else 0.0
            remaining = len(b_candidates) - idx
            eta = (remaining / b_rate) if b_rate > 1e-9 else float("inf")
            eta_str = f"{eta:.0f}s" if eta != float("inf") else "∞"
            print(f"[mitm] {idx:,}/{len(b_candidates):,} b | b_rate={b_rate:.2f}/s | rhs_rate={rhs_rate:,.0f}/s | ETA~{eta_str} | sols={len(sols)}")

    if progress:
        dt = time.perf_counter() - t0
        b_rate = (b_done / dt) if dt > 1e-9 else 0.0
        rhs_rate = (rhs_scanned_total / dt) if dt > 1e-9 else 0.0
        print(f"[mitm] done: b_checked={b_done:,} time={dt:.2f}s b_rate={b_rate:.2f}/s rhs_scanned={rhs_scanned_total:,} rhs_rate={rhs_rate:,.0f}/s sols={len(sols)}")

    return sols

# ============================================================
# Console Terminal Workflow
# ============================================================
def normalize_solution_list(sols: List[Tuple[Tuple[int, ...], int]]) -> set:
    return set(sols)

def run_validation_mode() -> None:
    print("\n=== VALIDATION MODE ===")
    print("Goal: brute-force vs MITM on SMALL bounds and assert they match.\n")

    k = int(input("Exponent k (e.g., 2, 4, 5): ").strip())
    m = int(input("LHS term count m (Euler-style use m=k-1): ").strip())
    N = int(input("Max value N (small! e.g. 20..60): ").strip())

    brute_est = estimate_nondec_tuple_count(N, m) * max(1, N - 1)
    print(f"\n[info] rough brute tuple count scale ~ {brute_est:,} (very rough)")
    if brute_est > 20_000_000:
        print("[warning] This brute run may be slow. Consider smaller N or smaller m.\n")

    split = input(f"MITM split r (blank for auto; default floor(m/2)={m//2}): ").strip()
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
        store_all=True,
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

    if m == 8 and split_r == 4 and HAS_CPP:
        cpp = cpp_m8_solutions(k, N, min_b=2, max_solutions=max_solutions, print_every=0)
        set_cpp = set(cpp)

        if set_cpp != mitm_set:
            only_cpp = set_cpp - mitm_set
            only_py  = mitm_set - set_cpp
            print("FAIL: C++ and Python MITM disagree.")
            print("Only in C++:", list(only_cpp)[:5])
            print("Only in Python:", list(only_py)[:5])
            return

        if set_cpp != brute_set:
            only_cpp = set_cpp - brute_set
            only_br  = brute_set - set_cpp
            print("FAIL: C++ and brute disagree.")
            print("Only in C++:", list(only_cpp)[:5])
            print("Only in brute:", list(only_br)[:5])
            return

        print("PASS: brute == Python MITM == C++ MITM on this parameter set.")

def run_research_mode() -> None:
    print("\n=== RESEARCH MODE ===")
    k = int(input("Exponent k (e.g., 9): ").strip())
    m = int(input("LHS term count m (e.g., 8): ").strip())

    min_b_str = input("Minimum b to start from (blank for 2): ").strip()
    min_b = int(min_b_str) if min_b_str else 2
    max_b = int(input("Maximum b (max_b): ").strip())

    max_val_str = input("Max LHS term value (blank for auto): ").strip()
    if max_val_str:
        max_val = int(max_val_str)
    else:
        # Auto defaults to keep MITM feasible. You can raise this if you have RAM.
        if m == 8:
            max_val = min(max_b, 250)
        elif m == 6:
            max_val = min(max_b, 1000)
        else:
            max_val = max_b

    split_str = input("Split r (blank for auto, forces 4 when m=8): ").strip()
    split_r = int(split_str) if split_str else None

    primes_str = input("Sieve primes (comma list) or blank for none (recommended for k>=7): ").strip()
    primes = [int(x) for x in primes_str.split(",") if x.strip()] if primes_str else None

    bucket_str = input("Bucket modulus (blank for auto; 4096 or 16384 recommended): ").strip()
    bucket_mod = int(bucket_str) if bucket_str else None

    if split_r is None:
        split_r = 4 if m == 8 else (2 if m == 4 else (m // 2))
    if bucket_mod is None:
        bucket_mod = 16384 if (m == 8 and split_r == 4 and (k >= 7 or max_b >= 200)) else (4096 if (k >= 7 or max_b >= 200) else None)

    store_all = input("Store ALL witnesses per sum? (y/N): ").strip().lower() == "y"
    grouping = input("Group scaled multiples into primitive classes? (y/N): ").strip().lower() == "y"
    max_s = input("Max solutions to collect (default 20): ").strip()
    max_solutions = int(max_s) if max_s else 20

    per_b_str = input("Print per-b line every how many b? (default 1; try 10 or 25 for speed): ").strip()
    per_b_every = int(per_b_str) if per_b_str else 1
        
    print("\n[note] This version uses LEGIT RHS cutoff (sorted by sum).")

    USE_CPP_THRESHOLD_K = 6
    use_cpp = HAS_CPP and (k >= USE_CPP_THRESHOLD_K)

    # For k>=8 and m in {8,9}, the classic 4+4 MITM needs an O(max_val^4) table and will blow up RAM.
    # This C++ backend is exact and RAM-safe: DFS + multi-prime pruning + 2-sum completion.
    if use_cpp and (k >= 8) and (m in (8, 9)):
        print(f"[backend] Using C++ RAM-safe DFS (m={m})")
        sieve_list = primes if primes is not None else []
        if m == 8:
            sols = cpp_m8_solutions(
                k=k,
                max_b=max_b,
                min_b=min_b,
                max_solutions=max_solutions,
                primitive_grouping=grouping,
                max_val=max_val,
                sieve_primes=sieve_list,
            )
        else:
            sols = cpp_m9_solutions(
                k=k,
                max_b=max_b,
                min_b=min_b,
                max_solutions=max_solutions,
                primitive_grouping=grouping,
                max_val=max_val,
                sieve_primes=sieve_list,
            )

    elif use_cpp and (m == 6 and split_r == 3):
        print("[backend] Using C++ RAM-safe DFS (m=6)")
        sieve_list = primes if primes is not None else []
        sols = cpp_m6_solutions(
            k=k,
            max_b=max_b,
            min_b=min_b,
            max_solutions=max_solutions,
            primitive_grouping=grouping,
            max_val=max_val,
            sieve_primes=sieve_list,
        )

    else:
        print("[backend] Using Python MITM")
        sols = mitm_solutions(
            k, m, max_b,
            min_b=min_b,
            split_r=split_r,
            primes=primes,
            bucket_mod=bucket_mod,
            store_all=store_all,
            max_solutions=max_solutions,
            progress=True,
            primitive_grouping=grouping,
            per_b_print_every=per_b_every,
        )


    print("\n=== OUTPUT ===")
    if not sols:
        print("No solutions found.")
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