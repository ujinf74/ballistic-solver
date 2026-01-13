<img width="2024" height="512" alt="banner" src="https://github.com/user-attachments/assets/f4e57e3f-f584-4938-a321-e9dd83dbbac3" />

**ballistic-solver** is a native C/C++ numerical solver that computes launch angles to intercept moving targets under gravity and quadratic air drag.
Unlike vacuum or closed-form solutions, this solver works in **real-time** and **fully nonlinear** conditions and converges even when trajectories are strongly curved.

---

## Quick start

### Python

```python
import ballistic_solver as bs

result = bs.solve(
    relPos0=(120, 30, 5),
    relVel=(2, -1, 0),
    v0=90,
    kDrag=0.002,
)

print(result["theta"], result["phi"], result["miss"])
```

### C ABI (single entry point)

```c
int32_t ballistic_solve(const BallisticInputs* in, BallisticOutputs* out);
```
See `ballistic_solver_c_api.h` for `BallisticInputs/Outputs` definitions.

---

## Demo (Unity)

Highly curved trajectories under strong air drag, still converging to a hit against moving targets.

https://github.com/user-attachments/assets/dcaf7479-cb94-477a-b71e-470a5b4c6004

---

## Why this solver

Many launch-angle solvers depend on vacuum assumptions or partially linearized models.
This project instead **simulates the projectile** and **solves the intercept numerically**, targeting robustness in real-time simulations and game/system integration.

---

## Key properties

* **Moving targets supported**
* **Strong air resistance supported**
* **Robust in strongly nonlinear regimes (no analytic assumptions)**
* **Best-effort result returned even without perfect convergence**
* **Explicit success / failure reporting**
* **Stable C ABI for multi-language use**
* **Header-only C++ core**
* **Low / High arc selection (v0.2.0)**

---

## Arc mode (v0.2.0)

The solver can be guided to converge to either the **low arc** or the **high arc** solution.

C ABI convention:

* `arcMode = 0` → Low
* `arcMode = 1` → High

### High arc example

https://github.com/user-attachments/assets/4334ed87-597e-4ad4-b21e-c1a1a17e8cd8

---

## C ABI: first-class interface

The core solver is written in C++, and exposed through a **stable C ABI** intended for long-term integration.

This enables usage from:

* C / C++
* Python (ctypes)
* C# / .NET / Unity (P/Invoke)
* Others via FFI

Prebuilt binaries are provided via GitHub Releases.

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

```
examples/dotnet/
```

On Windows, place `ballistic_solver.dll` next to the executable
(or ensure it is discoverable via PATH),
then call `ballistic_solve` via `DllImport`.

This works directly inside Unity.

---

## How it works (high level)

1. Simulate projectile motion using RK4 integration with drag
2. Track the closest approach between projectile and target
3. Express the miss at closest approach as an angular residual
4. Solve the nonlinear system using damped least squares (Levenberg–Marquardt)
5. Accelerate Jacobian updates with Broyden-style refinement
6. Return the best solution found

Failure cases are explicitly detected and reported.

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
