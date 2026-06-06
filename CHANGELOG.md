# Changelog

All notable changes to this project are documented here. Versions follow
`vMAJOR.MINOR.PATCH` git tags; the build derives its version from the tag.

## 0.7.0

### Breaking (C ABI v4)

- `BallisticInputs` gains an `int32_t preset` field (in the slot previously
  reserved as `_pad0`; struct layout is unchanged). `ABI` version bumped `3` -> `4`.
- `ballistic_inputs_apply_preset` now records the chosen preset, and
  `ballistic_solve` expands it through `make_params_preset`. Previously the
  C ABI carried only `dt`/`tMax`/`tolMiss`/`maxIter`, silently dropping the
  deeper solver tuning (line-search tries, lambda tries, finite-difference step,
  golden-section settings); presets are now fully applied across the C/.NET path.
- Callers must call `ballistic_inputs_init` (it sets `preset = 1`, Balanced);
  a zero-initialized struct now selects Fast (`preset = 0`).
- Added the `BALLISTIC_SOLVER_CALL` calling-convention macro to
  `ballistic_inputs_init`, `ballistic_accel_inputs_init`, and `ballistic_solve`
  declarations for consistency with the rest of the C ABI.

### Changed

- Solver: auxiliary multi-start fallback refactor (deduplicated residual path,
  clearer naming). High-arc and moving-target convergence verified at
  10000/10000 across low/high x stationary/moving case sets via
  `tools/bench_variants/compare_algorithms.py`.

### Documentation / tooling

- Documented the optional `relAcc` argument of the Python `solve()`.
- `benchmarks/linear_cases.py` now asserts a minimum success rate as a guard.
- Packaging: source distributions now include `LICENSE` and `SECURITY.md`.
