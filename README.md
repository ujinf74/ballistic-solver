# Official repository notice

This is the only official repository for ballistic-solver.

Official repository:
https://github.com/ujinf74/ballistic-solver

Do not download ZIP files, binaries, or installers from third-party copies, mirrors, or reuploads.
Use only the Releases page of this repository.

<img width="2024" alt="banner" src="https://github.com/user-attachments/assets/f4e57e3f-f584-4938-a321-e9dd83dbbac3" />

[![CI](https://github.com/ujinf74/ballistic-solver/actions/workflows/ci.yml/badge.svg?branch=main)](https://github.com/ujinf74/ballistic-solver/actions/workflows/ci.yml)
[![Release Native + PyPI](https://github.com/ujinf74/ballistic-solver/actions/workflows/release.yml/badge.svg)](https://github.com/ujinf74/ballistic-solver/actions/workflows/release.yml)
[![PyPI](https://img.shields.io/pypi/v/ballistic-solver)](https://pypi.org/project/ballistic-solver/)
[![License](https://img.shields.io/github/license/ujinf74/ballistic-solver)](https://github.com/ujinf74/ballistic-solver/blob/main/LICENSE)

**ballistic-solver** is a native C/C++ numerical solver that computes launch angles to intercept **moving targets** under **gravity** and **quadratic air drag**, with optional **wind**.

Unlike vacuum / closed-form solvers, this project **simulates the projectile** and **solves the intercept numerically**, aiming for robust real-time use even when trajectories are strongly curved.

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
solve(relPos0, relVel, v0, kDrag, arcMode=0, params=None) -> dict
```

* `relPos0`: target relative position at t=0 (x,y,z)
* `relVel`: target relative velocity (x,y,z)
* `v0`: muzzle speed (scalar)
* `kDrag`: quadratic drag coefficient
* `arcMode`: `0/1` or `"low"/"high"` (case-insensitive)
* `params`: optional `BallisticParams` for advanced tuning (gravity, wind, integrator and solver knobs)

Returned dict keys include:

* `success` (bool)
* `theta`, `phi` (radians)
* `miss` (closest-approach distance)
* `tStar` (time of closest approach)
* `relMissAtStar` (3-vector miss at `tStar`)
* `status` (SolveStatus integer)
* `message` (short diagnostic string)
* plus convergence diagnostics (`iterations`, `acceptedSteps`, `lastLambda`, `lastAlpha`)

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

1. Simulate projectile motion using RK4 integration with drag (+ wind)
2. Track the closest approach between projectile and target
3. Express the miss at closest approach as an angular residual
4. Solve the nonlinear system using damped least squares (Levenberg–Marquardt)
5. Accelerate Jacobian updates with Broyden-style refinement
6. Return the best solution found

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

Experimental solver-variant benchmarks live under `tools/bench_variants/`.
They are local development tools and are intentionally excluded from source distributions.

Benchmark numbers depend on CPU, OS, compiler, build type, Python version, and
whether native or Python entrypoints are measured. Record those fields when
publishing comparisons.

Reference result from a local Windows Release build, 500 generated linear-target
cases:

```text
fast:     median 0.094 ms, p95 0.228 ms, p95 miss 3.834e-02 m
balanced: median 0.182 ms, p95 0.452 ms, p95 miss 5.351e-03 m
precise:  median 0.199 ms, p95 0.569 ms, p95 miss 5.742e-06 m
```

High-arc moving-target convergence update, local Windows Release build, 500
generated cases:

| Solver configuration | Success | Median runtime | P95 runtime | P95 miss |
|---|---:|---:|---:|---:|
| Previous core path | 382/500 (76.40%) | 4.301 ms | 27.124 ms | 3.227e+02 m |
| v0.6.0 defaults | 490/500 (98.00%) | 1.845 ms | 2.535 ms | 8.809e-03 m |

---

## ABI notes

* Plain C layout across the ABI boundary
* Fixed-size arrays only
* No dynamic allocation across the boundary

---

## License

MIT License
