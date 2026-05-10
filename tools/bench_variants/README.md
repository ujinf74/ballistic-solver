# Solver Variant Benchmarks

This directory is for local solver experiments that are not part of the
distributed package.

## Scope

- `compare_algorithms.py` compares the released solver against experimental
  variants.
- `native_variants/` builds benchmark-only pybind11 variants such as `no_aux`
  and `aux_fixed_point`.
- `build_bench_variants.ps1` builds the native extension into the repo-local
  `.bench_resume_site/` import site.
- `compare_algorithms.py --target-mode stationary` benchmarks stationary targets.
- `compare_algorithms.py --min-hit-z 0` filters generated cases to nonnegative
  hit altitude.

The distributable benchmark surface stays in `benchmarks/linear_cases.py`.
