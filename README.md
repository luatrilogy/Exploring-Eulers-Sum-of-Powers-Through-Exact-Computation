# Exploring Euler’s Sum of Powers Through Exact Computation  
### Modular Arithmetic, Sumsets, and Search Design

This project is an **exact computational framework** for exploring Euler-type sum-of-powers Diophantine equations

\[
x_1^k + x_2^k + \cdots + x_m^k = b^k,\quad x_i \le b,
\]

motivated by **Euler’s Sum of Powers conjecture** (the boundary regime is typically \(m<k\), e.g. \(m=k-1\)). The code focuses on **mathematical legitimacy**: it avoids floating-point screening, uses exact integer arithmetic for verification, and only restricts the **b-range** (not the \(x_i\) domain) when requested.

In addition to boundary-style searches, I use even-split regimes such as \(m=6\) (split \(3+3\)) and \(m=8\) (split \(4+4\)) as **stress tests** for meet-in-the-middle performance and memory behavior.

---

## Key ideas

### 1) Meet-in-the-middle (MITM)
The core method precomputes sums of \(r\) \(k\)-th powers and searches against complementary \((m-r)\)-term sums to match \(b^k\). Nondecreasing tuples are used to reduce duplicates while preserving completeness.

### 2) Modular sieve on \(b\) using residue sumsets
For primes \(p\), compute the set of \(k\)-th power residues \(\{x^k \bmod p\}\), then compute the reachable residues formed by the **m-fold sumset**. This is implemented efficiently via **bitmask cyclic shifts**, allowing the search to discard any \(b\) whose \(b^k \bmod p\) is unreachable as a sum of \(m\) residues.

### 3) Exact verification + (optional) C++ backend
All candidate matches are verified with exact arithmetic. An optional C++/pybind11 extension (`euler_ext.cpp`) implements specialized MITM kernels for hard regimes (e.g., \(m=8\) split \(4+4\), \(m=6\) split \(3+3\)) and supports residue-indexed bucketing for speed.

---

## Repository structure
- `app.py` — Python utilities / helpers
- `euler_ext.cpp` — optional pybind11 C++ backend (exact MITM kernels)
- `setup.py` — build script for the C++ extension (requires editing include paths on Windows)
---

## Quickstart (Python-only)

**Requirements:** Python 3.8+ (recommended 3.10+)