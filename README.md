<img width="2024" alt="banner" src="https://github.com/user-attachments/assets/f4e57e3f-f584-4938-a321-e9dd83dbbac3" />

[![CI](https://github.com/ujinf74/ballistic-solver/actions/workflows/ci.yml/badge.svg?branch=main)](https://github.com/ujinf74/ballistic-solver/actions/workflows/ci.yml)
[![Release Native + PyPI](https://github.com/ujinf74/ballistic-solver/actions/workflows/release.yml/badge.svg)](https://github.com/ujinf74/ballistic-solver/actions/workflows/release.yml)
[![PyPI](https://img.shields.io/pypi/v/ballistic-solver)](https://pypi.org/project/ballistic-solver/)
[![License](https://img.shields.io/github/license/ujinf74/ballistic-solver)](https://github.com/ujinf74/ballistic-solver/blob/main/LICENSE)

**ballistic-solver** is a native C/C++ numerical solver that computes the launch angles to hit **moving targets** under **gravity** and **quadratic air drag**, with optional **wind**. Unlike vacuum / closed-form solvers, it **simulates the projectile** and **solves the intercept numerically**, so it stays accurate even when trajectories are strongly curved. Ships as a header-only C++ core with a stable C ABI and Python (`pip install ballistic-solver`), C#/Unity, and Godot bindings.

https://github.com/user-attachments/assets/7de04137-61d9-4cf0-be5e-9804d6a9c67b

> Real-time 3D showcase — the built-in tracker leading a noisy, diving target with range rings, tracer trails, and a live hit-rate HUD (`examples/godot/`).

**Applications.** The solver only sees relative kinematics, so the same code
drops into many domains — most visibly **games and simulation** (turret
lead-targeting, projectile interception; see the demo above), and equally
**robotics** (interception / capture), tracking, and other aiming / targeting
problems. An optional position-only tracker (`bs.TargetTracker`) adds lead-fire
from a noisy track. See `benchmarks/intercept_eval.py` and `examples/viz/` for
cross-domain examples and accuracy comparisons.

---

## Quick start

### Python (PyPI)

```bash
pip install ballistic-solver
```

Requires Python **>= 3.10**.

```python
import ballistic_solver as bs

result = bs.solve(
    relPos0=(120, 30, 5),
    relVel=(2, -1, 0),
    v0=90,
    kDrag=0.002,
)

print(result["theta"], result["phi"], result["miss"])
print(result["success"], result["status"], result["message"])
```

For tighter convergence without manually tuning every knob:

```python
params = bs.params_preset("precise")
result = bs.solve((120, 30, 5), (2, -1, 0), 90, 0.002, params=params)
```

For repeated solves, Python also provides a thin convenience wrapper:

```python
solver = bs.Solver.preset("precise")
result = solver.solve((120, 30, 5), (2, -1, 0), 90, 0.002)
```

---

## Demo (Unity)

Highly curved trajectories under strong air drag, still converging to a hit against moving targets.

https://github.com/user-attachments/assets/c0c69cdd-0dd4-4606-9c7d-f21dd002d7f7

---

## Why this solver

Many launch-angle solvers depend on vacuum assumptions or partially linearized models.
This project instead **simulates the projectile** and **solves the intercept numerically**, targeting robustness in real-time simulations and integration scenarios.

---

## Key properties

* Moving targets supported
* Constant-acceleration targets supported via the extended API
* Strong air resistance (quadratic drag) supported
* Low / High arc selection (since v0.2)
* Wind vector supported (since v0.3)
* Extended C ABI utilities (since v0.4)
* Quartic vacuum-lead initialization for moving targets (since v0.6)
* High-arc moving-target auxiliary multi-start fallback (since v0.6.1)
* Coordinate-residual Gauss–Newton core with vacuum-seeded Jacobian (default since v1.0; ~15–25% faster, auxiliary-residual method kept as `solve_aux`)
* Fast / balanced / precise solver presets
* Physical drag helper: `kDrag = 0.5 * rho * Cd * area / mass`
* Robust in strongly nonlinear regimes (no analytic assumptions)
* Best-effort result returned even without perfect convergence
* Explicit success / failure reporting (+ diagnostic message)
* Stable C ABI for multi-language use
* Header-only C++ core
* Easy install via PyPI: `pip install ballistic-solver`

---

## Python API

### `solve(...)`

```python
solve(relPos0, relVel, v0, kDrag, arcMode=None, params=None, relAcc=None) -> dict
```

* `relPos0`: target relative position at t=0 (x,y,z)
* `relVel`: target relative velocity (x,y,z)
* `v0`: muzzle speed (scalar)
* `kDrag`: quadratic drag coefficient
* `arcMode`: `0/1` or `"low"/"high"` (case-insensitive); `None` keeps the `params` value
* `params`: optional `BallisticParams` for advanced tuning (gravity, wind, integrator and solver knobs)
* `relAcc`: optional constant target relative acceleration (x,y,z); same as `solve_accel`

Returned dict keys include:

* `success` (bool)
* `theta`, `phi` (radians)
* `miss` (closest-approach distance)
* `tStar` (time of closest approach)
* `relMissAtStar` (3-vector miss at `tStar`)
* `status` (SolveStatus integer)
* `message` (short diagnostic string)
* plus convergence diagnostics (`iterations`, `acceptedSteps`, `lastLambda`, `lastAlpha`)

`solve` uses the coordinate-residual core (Gauss–Newton on the 3D closest-approach
miss, vacuum-seeded Jacobian, multistart fallback). It is the same signature as
before and a drop-in replacement, typically ~15–25% faster at equal robustness.

### `solve_aux(...)`

```python
solve_aux(relPos0, relVel, v0, kDrag, arcMode=None, params=None, relAcc=None) -> dict
```

The previous default: the auxiliary-residual method (re-aim through the vacuum
inverse, solved by damped least squares with a multistart fallback). Same return
dict as `solve`. Kept for compatibility and for reproducing the auxiliary-residual
results.

### `solve_predicted(...)`

```python
solve_predicted(predictor, v0, kDrag, arcMode=None, params=None) -> dict
```

Lead-fire against a position-only trajectory predictor instead of explicit
`relVel`/`relAcc`. `predictor(t)` returns the predicted relative position `(x, y, z)`
at time `t >= 0` (e.g. from a Kalman/IMM track). Internally a local
constant-acceleration model is fitted at the intercept time and iterated to a
fixed point — exact for constant-velocity / constant-acceleration predictors and
second-order accurate at the intercept for smoothly curving tracks. Returns the
same dict as `solve`.

### `TargetTracker`

```python
tracker = bs.TargetTracker(processNoise=1.0, measNoise=0.25)
for t, pos in measurements:        # (timestamp, relative position)
    tracker.update(t, pos)
lead = tracker.predict(1.5)        # predicted relative position 1.5 s ahead
result = tracker.solve(v0=95.0, kDrag=0.0015)   # lead-fire via the predictor seam
```

A position-only constant-acceleration Kalman tracker (per-axis, white-noise-jerk
model). It denoises the track and recovers velocity/acceleration without
finite-difference noise, then feeds the predictor seam — so callers never supply
`relVel`/`relAcc`. See `benchmarks/predictor_eval.py` for a lead-prediction
accuracy comparison against a constant-velocity baseline.

The real-time 3D showcase at the top of this README is built on this tracker —
the full project is in `examples/godot/` (default scene `demo/intercept_demo.tscn`).

### Utilities

```python
params = bs.params_preset("fast")  # or "balanced", "precise"
kDrag = bs.k_drag_from_physical(
    airDensity=1.225,
    dragCoefficient=0.30,
    area=0.00426,
    mass=0.145,
)
```

For constant-acceleration targets:

```python
result = bs.solve_accel(
    relPos0=(120, 30, 5),
    relVel=(2, -1, 0),
    relAcc=(0, 0.2, 0),
    v0=90,
    kDrag=0.002,
    params=bs.params_preset("precise"),
)
```

---

## Solver validity note

The solver internally integrates projectile motion using:

* 4th-order Runge–Kutta (RK4)
* Fixed timestep `dt`
* Quadratic drag: `a = (0, 0, -g) - kDrag * |v - wind| * (v - wind)`
* Wind as air velocity

To match in-game ballistics, your runtime simulation must use the **same physical model and integrator configuration**.

If your game uses a different integrator (e.g., Euler) or a different timestep, the computed launch angles may not hit even if the solver reports success.

---

## C ABI (stable interface)

Primary intercept API:

```c
void ballistic_inputs_init(BallisticInputs* in);
int32_t ballistic_solve(const BallisticInputs* in, BallisticOutputs* out);
```

Return value policy:

* `0`: API call completed and `out` was filled.
* `<0`: API-level failure, such as null pointers or an internal exception.
* Numerical solve success is reported separately by `out->success` and `out->status`.
* ABI v3 adds convergence diagnostics to `BallisticOutputs`: `iterations`,
  `acceptedSteps`, `lastLambda`, and `lastAlpha`.
* ABI v4 adds a `preset` field to `BallisticInputs` (0=Fast, 1=Balanced,
  2=Precise). `ballistic_inputs_apply_preset` records it so that the full preset
  tuning (line-search, lambda tries, finite-difference step, golden-section)
  is carried into `ballistic_solve`, not only `dt`/`tMax`/`tolMiss`/`maxIter`.

Callers should check both:

```c
int32_t rc = ballistic_solve(&in, &out);
if (rc != 0) {
    /* API call failed */
}
if (!out.success) {
    /* Solver ran but did not satisfy the requested tolerance. */
}
```

Since **v0.4.0**, additional utility functions are available:

```c
void ballistic_rk4_step(...);

int32_t ballistic_simulate_trajectory(...);
int32_t ballistic_simulate_trajectory_from_angles(...);

int32_t ballistic_find_closest_approach(...);

int32_t ballistic_vacuum_arc_angles_to_point(...);
void ballistic_initial_guess_vacuum_lead(...);
```

Additional extended APIs include:

```c
void ballistic_accel_inputs_init(BallisticAccelInputs* in);
int32_t ballistic_solve_accel(const BallisticAccelInputs* in, BallisticOutputs* out);
int32_t ballistic_inputs_apply_preset(BallisticInputs* in, int32_t preset);
int32_t ballistic_k_drag_from_physical(...);
int32_t ballistic_make_relative_motion(...);
```

See `ballistic_solver_c_api.h` for full signatures and parameter definitions.

This enables usage from:

* C / C++
* Python (ctypes via the C ABI)
* C# / .NET / Unity (P/Invoke)
* Unity (UnityEngine.Vector3 wrapper example)
* Godot 4 (GDExtension addon-style example)
* Others via FFI

Prebuilt native binaries are provided via GitHub Releases.

---

## Arc mode (since v0.2)

C ABI convention:

* `arcMode = 0` → Low
* `arcMode = 1` → High

High arc example:

https://github.com/user-attachments/assets/4334ed87-597e-4ad4-b21e-c1a1a17e8cd8

---

## Wind (since v0.3)

C ABI convention:

* `wind[3]` = air velocity vector (same frame as `relPos0/relVel`)
* Drag uses relative airspeed: `v_rel = v_projectile - wind`

Wind demo:

https://github.com/user-attachments/assets/1cd998cf-34db-4a74-8817-c6393227ef4e

---

## Using prebuilt binaries (C ABI)

Download the archive for your platform from **Releases**.

Each release contains:

* Shared library

  * Windows: `ballistic_solver.dll`
  * Linux: `libballistic_solver.so`
  * macOS: `libballistic_solver.dylib`
* C ABI header: `ballistic_solver_c_api.h`

---

## C# / Unity usage

A C# P/Invoke example is available in:

```text
examples/dotnet/
```

On Windows, place `ballistic_solver.dll` next to the executable
(or ensure it is discoverable via PATH),
then call `ballistic_solve` via `DllImport`.

This works directly inside Unity.

---

## How it works (high level)

1. Build a vacuum-lead initial guess from the moving-target quartic.
2. Seed the inverse Jacobian analytically from the vacuum-arc map (no trajectory
   integration), as a cheap preconditioner.
3. Simulate projectile motion using RK4 integration with drag (+ wind) and track
   the closest approach between projectile and target.
4. Use the 3D closest-approach miss vector itself as the residual (it is
   well-conditioned), and solve the launch angles by Gauss-Newton on it.
5. Refine the Jacobian with Broyden-style rank-1 updates; take full steps
   (no Levenberg–Marquardt damping, no line search).
6. If the warm start does not reach tolerance, retry from a small arc-appropriate
   theta-grid multistart.
7. Return the best solution found.

`solve` uses this coordinate-residual core. The earlier auxiliary-residual method
(re-aim through the vacuum inverse, with damped least squares) remains available
as `solve_aux` for compatibility and reproducibility.

Failure cases are explicitly detected and reported.

---

## When solving can fail

The solver returns a best-effort result even when it cannot satisfy `tolMiss`.
Common causes include:

* invalid inputs (`v0 <= 0`, non-positive `g`, `dt`, `tMax`, or `maxIter`)
* geometrically unreachable vacuum targets
* targets that require an intercept beyond `tMax`
* strong drag or high-arc cases where the selected arc branch cannot be maintained
* iteration limits that are too tight for the requested tolerance

Use `success`, `status`, `message`, and `miss` together when deciding whether to
accept a solution.

---

## Status codes (SolveStatus)

`BallisticOutputs.status` / Python `result["status"]` corresponds to:

* `0` = Ok
* `1` = InvalidInput
* `2` = InitialResidualFailed
* `3` = JacobianFailed
* `4` = LMStepSingular
* `5` = ResidualFailedDuringSearch
* `6` = LineSearchRejected
* `7` = LambdaTriesExhausted
* `8` = MaxIterReached

The default coordinate-residual `solve` emits `Ok`, `InvalidInput`, `JacobianFailed`,
or `MaxIterReached`. The Levenberg–Marquardt codes (`4`–`7`) are specific to the
auxiliary-residual `solve_aux` path.

`message` contains a short diagnostic string.

---

## Build from source

```bash
cmake -S . -B build
cmake --build build -j
ctest --test-dir build
```

The shared library target is `ballistic_solver`.

## Regression and benchmark

Distributed regression and benchmark scripts are available:

```bash
python tests/random_regression.py
python benchmarks/linear_cases.py
```

The CI workflow runs CTest smoke coverage plus the Python regression script.
The regression script includes analytic vacuum checks, constructed moving-target
vacuum cases, unreachable-target checks, randomized linear cases, and
constant-acceleration API smoke coverage.

Benchmark numbers depend on CPU, OS, compiler, build type, Python version, and
whether native or Python entrypoints are measured. The local reference below was
measured on Windows NT 10.0.26200.0, AMD Ryzen 9 6900HS, Python 3.12.8, MSVC
14.44.35207, x64 Release build, through the Python native extension.

Reference preset result from the same local Release build:

```text
fast:     500/500, median 0.054 ms, p95 0.237 ms, p95 miss 2.810e-02 m
balanced: 500/500, median 0.115 ms, p95 0.492 ms, p95 miss 7.053e-03 m
precise:  500/500, median 0.168 ms, p95 0.669 ms, p95 miss 5.596e-06 m
```

Current default solver (coordinate-residual core) results from 10,000 generated
cases with seed `20260503`. The last row is the high-arc moving-target stress case:

| Case set | Success | Median runtime | P95 runtime | P95 miss |
|---|---:|---:|---:|---:|
| Low arc, moving target | 10000/10000 (100.00%) | 0.168 ms | 0.510 ms | 6.450e-03 m |
| High arc, stationary target | 10000/10000 (100.00%) | 0.785 ms | 1.604 ms | 8.018e-03 m |
| High arc, moving target | 10000/10000 (100.00%) | 1.126 ms | 2.350 ms | 8.341e-03 m |

---

## ABI notes

* Plain C layout across the ABI boundary
* Fixed-size arrays only
* No dynamic allocation across the boundary

---

## Official repository

The only official repository is https://github.com/ujinf74/ballistic-solver.
Install from PyPI (`pip install ballistic-solver`) or this repository's Releases
page — not from third-party copies, mirrors, or reuploads.

---

## License

MIT License
