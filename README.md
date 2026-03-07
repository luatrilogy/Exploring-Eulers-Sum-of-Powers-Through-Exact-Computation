# Exploring Euler's Sum of Powers Through Exact Computation
### Modular Arithmetic, Sumsets, and Search Design

![Python](https://img.shields.io/badge/python-3.8%2B-blue.svg)
![C++ Backend](https://img.shields.io/badge/C%2B%2B-backend%20optional-00599C.svg)
![Status](https://img.shields.io/badge/status-research%20prototype-orange.svg)

Abstract: https://drive.google.com/file/d/1_poJwWQcWq1QKoYczQnibSPVc94Suikj/view?usp=sharing

This project is an **exact computational framework** for exploring Euler-type sum-of-powers Diophantine equations motivated by **Euler's Sum of Powers conjecture**. The code focuses on **mathematical legitimacy**: it avoids floating-point screening, uses exact integer arithmetic for verification, and only restricts the **b-range** when requested.

In addition to boundary-style searches, I use even-split regimes such as `(m = 6)` split `(3+3)` and `(m = 8)` split `(4+4)` as **stress tests** for meet-in-the-middle performance and memory behavior.

---

## Table of contents
- [Key ideas](#key-ideas)
- [Repository structure](#repository-structure)
- [Quickstart (Python-only)](#quickstart-python-only)
- [Optional C++ backend](#optional-c-backend)
- [Troubleshooting](#troubleshooting)

---

## Key ideas

### 1) Meet-in-the-middle (MITM)
The core method precomputes sums of `r` `k`-th powers and searches against complementary `(m-r)`-term sums to match `b^k`. Nondecreasing tuples are used to reduce duplicates while preserving completeness.

### 2) Modular sieve on `b` using residue sumsets
For primes `p`, compute the set of `k`-th power residues `{x^k mod p}`, then compute the reachable residues formed by the **m-fold sumset**. This is implemented efficiently via **bitmask cyclic shifts**, allowing the search to discard any `b` whose `b^k mod p` is unreachable as a sum of `m` residues.

### 3) Exact verification + (optional) C++ backend
All candidate matches are verified with exact arithmetic. An optional C++/pybind11 extension (`euler_ext.cpp`) implements specialized MITM kernels for hard regimes (for example `(m = 8)` split `(4+4)` and `(m = 6)` split `(3+3)`), plus DFS + 2-sum routines for larger-memory regimes.

---

## Repository structure
- `app.py` - interactive solver entry point (validation mode + research mode)
- `euler_ext.cpp` - optional pybind11 C++ backend (exact MITM kernels)
- `setups.py` - build script for the C++ extension (edit include paths on Windows as needed)
- `visualization.py` - Manim scenes used to present the problem, method, and systems insights
- `media/` - generated Manim outputs (videos, SVG/TEX artifacts)
- `old versions/` - archived prototypes and experiments

---

## Quickstart (Python-only)

**Requirements:** Python 3.8+ (recommended 3.10+)

Create and activate a virtual environment (recommended):

```powershell
python -m venv .venv
.\.venv\Scripts\Activate.ps1
```

Install packages:

```powershell
python -m pip install --upgrade pip
python -m pip install tqdm pybind11 manim
```

Run the Python solver:

```powershell
python app.py
```

The menu includes:
- Validation mode: brute force vs Python MITM (and C++ comparison when available)
- Research mode: long-run exact searches with optional sieve and bucketing
- Quick known identity check for `27^5 + 84^5 + 110^5 + 133^5 = 144^5`

---

## Optional C++ backend

If you want the accelerated backend, first verify include paths in `setups.py` (especially pybind11/Boost/vcpkg paths for your machine), then build:

```powershell
python setups.py build_ext --inplace
```

If the extension loads, `app.py` will automatically route eligible workloads to C++.

Quick import check:

```powershell
python -c "import euler_ext; print('C++ extension loaded')"
```

---

## Troubleshooting
- `ImportError: No module named euler_ext`: run `python setups.py build_ext --inplace` and confirm build paths in `setups.py`.
- Build errors on Windows: check compiler toolchain, pybind11 install, and vcpkg include path configuration.
- Manim render fails: ensure a full Manim installation (including LaTeX/FFmpeg dependencies required by your setup).
- Slow runs in research mode: use sieve primes, reduce `max_val`, or switch to C++ backend for supported `(m, r)` splits.