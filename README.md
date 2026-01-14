<img width="2024" height="512" alt="banner" src="https://github.com/user-attachments/assets/f4e57e3f-f584-4938-a321-e9dd83dbbac3" />

**ballistic-solver** is a native C/C++ numerical solver that computes launch angles to intercept **moving targets** under **gravity** and **quadratic air drag**.
Unlike vacuum / closed-form approaches, this solver **simulates the projectile** and **solves the intercept numerically**, targeting robustness for real-time simulations even when trajectories are strongly curved.

---

## Quick start

### Python (wheel)

```bash
pip install ballistic-solver
````

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

#### Signature (actual binding)

```python
solve(relPos0, relVel, v0, kDrag, arcMode=0, params=None) -> dict
```

* `arcMode` accepts `0/1` or `"low"/"high"` (case-insensitive).
* **Important:** the `arcMode` argument **overrides** `params.arcMode` even if `params` is provided.

Returned dict keys include:

* `success`, `theta`, `phi`, `miss`, `tStar`, `relMissAtStar`
* `status`, `message`
* `iterations`, `acceptedSteps`, `lastLambda`, `lastAlpha`

### C ABI (stable interface)

```c
void ballistic_inputs_init(BallisticInputs* in);
int32_t ballistic_solve(const BallisticInputs* in, BallisticOutputs* out);
```

See `ballistic_solver_c_api.h` for `BallisticInputs/Outputs` definitions and defaults.

---

## Demo (Unity)

Highly curved trajectories under strong air drag, still converging to a hit against moving targets.

[https://github.com/user-attachments/assets/dcaf7479-cb94-477a-b71e-470a5b4c6004](https://github.com/user-attachments/assets/dcaf7479-cb94-477a-b71e-470a5b4c6004)

---

## Key properties

* Moving targets supported
* Strong air resistance (quadratic drag) supported
* Robust in strongly nonlinear regimes (no analytic assumptions)
* Best-effort result returned even without perfect convergence
* Explicit success / failure reporting (+ diagnostic message)
* Stable C ABI for multi-language use
* Header-only C++ core
* Low / High arc selection (since v0.2.0)

---

## Arc mode (since v0.2.0)

C ABI convention:

* `arcMode = 0` → Low
* `arcMode = 1` → High

High arc example:

[https://github.com/user-attachments/assets/4334ed87-597e-4ad4-b21e-c1a1a17e8cd8](https://github.com/user-attachments/assets/4334ed87-597e-4ad4-b21e-c1a1a17e8cd8)

---

## How it works (high level)

1. Simulate projectile motion using RK4 integration with drag (+ wind in core/C ABI)
2. Track the closest approach between projectile and target
3. Express the miss at closest approach as an angular residual
4. Solve the nonlinear system using damped least squares (Levenberg–Marquardt)
5. Accelerate Jacobian updates with Broyden-style refinement
6. Return the best solution found

Failure cases are explicitly detected and reported.

---

## Advanced parameters (C ABI)

`BallisticInputs` exposes additional controls (all have defaults via `ballistic_inputs_init`):

* `g` (gravity)
* `wind[3]` (wind vector)
* `dt` (RK4 step)
* `tMax` (max sim time)
* `tolMiss` (hit tolerance)
* `maxIter` (LM iterations)

Note: The current Python binding exposes many tuning parameters via `BallisticParams`,
but does **not** expose `wind` as a Python field.

---

## Status codes (SolveStatus)

`BallisticOutputs.status` / Python `result["status"]` corresponds to the core enum:

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

## Using prebuilt binaries

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

---

## Build from source

```bash
cmake -S . -B build
cmake --build build -j
ctest --test-dir build
```

The shared library target is `ballistic_solver`.

---

## ABI notes

* Plain C layout across the ABI boundary
* Fixed-size arrays only
* No dynamic allocation across the boundary

---

## License

MIT License
