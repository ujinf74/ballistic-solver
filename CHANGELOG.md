# Changelog

All notable changes to this project are documented here. Versions follow
`vMAJOR.MINOR.PATCH` git tags; the build derives its version from the tag.

## 1.2.0

### Added

- Wheels for CPython 3.13, 3.14 and 3.15, and for Linux `aarch64`, musllinux
  (Alpine) and Windows `ARM64`. Only cp310-cp312 on x86_64 Linux, AMD64
  Windows and macOS universal2 were published before, so `pip install` on 3.13
  or newer fell back to an sdist build and failed for anyone without CMake and
  a C++17 compiler. The ARM wheels are built on native ARM runners rather than
  under emulation, so their test command runs as well.

### Changed

- cibuildwheel 3.3.1 -> 4.3.2, and the pybind11 floor raised to `>=3.0`, which
  CPython 3.15 requires. `requires-python` stays at `>=3.10`.

### Fixed

- Release wheel artifacts were named after `runner.os`, which collides as soon
  as two runners of the same OS build different architectures;
  `upload-artifact@v4` rejects the duplicate name. Keyed on `matrix.os`
  instead.

### Documentation / tooling

- README corrections: the FetchContent example was pinned to the superseded
  v1.0.0; the Python section gave only a `>= 3.10` floor rather than the wheel
  matrix; the native archives are single-architecture builds whose filenames
  do not say so, so someone on an Intel Mac downloads an arm64-only dylib; and
  `lastLambda` / `lastAlpha` / `acceptedSteps` were presented as convergence
  diagnostics without noting that the default coordinate-residual core never
  enters a Levenberg-Marquardt path and leaves them at `NaN`, `NaN` and `0`.
- The `TargetTracker` snippet looped over an undefined `measurements`, so it
  could not be pasted and run.
- Trimmed the Key properties list from 17 bullets to 9, dropping a version
  chronology that duplicated this changelog, internals documented elsewhere,
  and bullets restating the comparison table.

## 1.1.0

### Added

- Modern C++ API surface in `<ballistic_solver.hpp>`: `bs::solve` / `bs::solve_aux`
  taking `bs::Problem` and `bs::Options` and returning `bs::Intercept`, namespaced
  and compiled into the shared library so no internal symbols leak into consumers.
  The C ABI (`ballistic_solve`) remains the stable / FFI boundary.
- Test coverage for that C++ surface (`tests/smoke_test_cpp.cpp`, registered with
  CTest): `solve` / `solve_aux`, the preset and arc options, the `rel_acc` path
  and the rejected-input path. `examples/cpp/basic.cpp` is now a CMake target as
  well, so neither can rot unnoticed.

### Changed

- Decomposed the two oversized solver functions: `solve_launch_angles`
  (328 -> 220 lines) and `solve_auxiliary_multistart` (232 -> 44), hoisting the
  finite-difference Jacobian that had been duplicated between them into a single
  `jacobian_angles_fd`, and splitting the multistart into
  `collect_multistart_seeds` and `run_multistart_seed`. Pure extraction: no
  operation was reordered and no public signature changed. Verified bit-for-bit
  identical output over a deterministic 9000-solve sweep against the
  pre-refactor baseline.
- `examples/cpp/basic.cpp` now uses the public `bs::` API instead of reaching
  into the internal core headers, matching what the README documents.

### Fixed

- The C smoke test was compiled away under Release (`NDEBUG`), silently
  disabling every assertion. Worse, five API calls sat inside the assert
  expressions, so `ballistic_k_drag_from_physical`,
  `ballistic_inputs_apply_preset`, the precise-preset solve and
  `ballistic_solve_accel` were never invoked in CI at all. Checks are now
  explicit and survive `NDEBUG`.
- That also surfaced a stale assertion on `BallisticOutputs.lastLambda`.
  `lastLambda` / `lastAlpha` are Levenberg-Marquardt diagnostics, and the
  default core (coordinate-residual) uses neither damping nor a line search, so
  it leaves both at NaN by design. Consumers should not read them on the default
  path.

### Documentation / tooling

- PyPI package metadata: author, keywords, trove classifiers and project URLs.
  The published package previously had none of these, so it was unreachable
  through PyPI search and its page did not link back to the repository.
- `SECURITY.md` now names 1.0.x as the supported release line (it still listed
  0.7.x).
- Added `.mailmap` so `git log` / `git shortlog` report a single author identity
  across the history.

## 1.0.0

### Added

- Target predictor seam: `solve_launch_angles_predicted(TargetPredictor, ...)`
  (C++, `include/bs/predictor.hpp`) and `solve_predicted(predictor, v0, kDrag, ...)`
  (Python) accept a position-only trajectory predictor `f(t)` instead of explicit
  `relVel`/`relAcc`. Implemented additively over the existing solver by fitting a
  local constant-acceleration model at the intercept time and iterating it to a
  fixed point: exact for constant-velocity / constant-acceleration predictors and
  second-order accurate at the intercept for smoothly curving tracks. The core
  solver is unchanged.
- `TargetTracker` (C++ `include/bs/kalman.hpp`, Python `bs.TargetTracker`): a
  position-only constant-acceleration Kalman tracker (per-axis, white-noise-jerk
  model). Ingests `(timestamp, relative position)` and exposes `predict(tau)` and
  a one-call `solve(v0, kDrag, ...)` via the predictor seam, so callers need not
  supply `relVel`/`relAcc`. On a 0.5 m-noise, 20 Hz feed it reduces 1 s lead-
  prediction RMS by ~4-8x versus a two-point constant-velocity baseline on
  constant-velocity, accelerating, and coordinated-turn tracks
  (`benchmarks/predictor_eval.py`).

- Examples: the Godot GDExtension exposes a `BallisticTracker` class (wrapping
  the C++ `TargetTracker` directly), and a new `demo/intercept_demo.tscn`
  showcase tracks a moving target from a noisy position-only feed and fires a
  scored burst. `benchmarks/intercept_eval.py` covers robotics / security /
  agriculture / defense scenarios, and `examples/viz/intercept_viz.py` renders
  miss-vs-range curves and an engagement GIF (matplotlib, optional).
- Coordinate-residual solver core (`include/bs/solve_coord_lead.hpp`): solves the
  launch angles by Gauss-Newton on the 3D closest-approach miss vector, with the
  inverse Jacobian seeded analytically from the vacuum-arc map (no trajectory
  integration) and refined by Broyden, plus an arc-grid multistart fallback. No
  Levenberg-Marquardt damping or line search.
- `solve_aux` (Python): the prior auxiliary-residual solver, kept for
  compatibility and for reproducing the auxiliary-residual results.

### Changed

- Default solver swapped to the coordinate-residual core: `solve` (Python) and
  `ballistic_solve` / `ballistic_solve_accel` (C ABI) now use it. Same signatures
  and result dict; ~15-25% faster at equal (100%) success across the low/high x
  moving/stationary case sets, and simpler (drops the auxiliary residual, LM, and
  line search from the default path). The auxiliary-residual method remains as
  `solve_aux`; `SolveStatus` codes 4-7 are specific to it. C ABI struct layout is
  unchanged (binary-compatible).
- Godot `intercept_demo`: round tracer trails, lead-to-intercept line, a denoised
  track-estimate marker, radar range rings, and a 3/4 camera; the
  `BallisticTracker` firing solve uses the precise preset.
- Refactor: split `ballistic_solver_core.hpp` into focused headers under
  `include/bs/` (vec3, math_utils, params, integrator, target, closest_approach,
  vacuum_lead, residual, lm, solve, predictor). Behavior-preserving; the umbrella
  header keeps the existing include surface.

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
  10000/10000 across low/high x stationary/moving case sets (local
  algorithm-comparison harness; `tests/random_regression.py` guards a subset).

### Documentation / tooling

- Documented the optional `relAcc` argument of the Python `solve()`.
- `benchmarks/linear_cases.py` now asserts a minimum success rate as a guard.
- Packaging: source distributions now include `LICENSE` and `SECURITY.md`.
