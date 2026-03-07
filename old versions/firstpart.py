import time
import math
from typing import Dict, List, Tuple, Iterable, Optional, Union

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

    max_s = input("Max solutions to collect (default 50): ").strip()
    max_solutions = int(max_s) if max_s else 50

    print("\nRunning brute force...")
    brute = brute_solutions(k, m, N, max_solutions=max_solutions)

    if brute:
        print("\nExample solutions (lhs -> b):")
        for lhs, b in brute[:min(10, len(brute))]:
            print(f"  {lhs} -> {b}")


def main():
    while True:
        print("\n==============================")
        print(" Euler / Sum-of-Powers Search ")
        print("==============================")
        print("1) Validation mode (brute vs MITM)")
        #print("2) Research mode (MITM + sieve/buckets + packed tables)")
        #print("3) Quick known check: 27^5+84^5+110^5+133^5 == 144^5")
        #print("4) Exit")
        choice = input("Select: ").strip()

        if choice == "1":
            run_validation_mode()
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

    # ================================
# CONSOLE UI: Validation + Research Modes
# ================================