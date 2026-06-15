# Changelog

All notable changes to this project are documented here. Versions follow
`vMAJOR.MINOR.PATCH` git tags; the build derives its version from the tag.

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
