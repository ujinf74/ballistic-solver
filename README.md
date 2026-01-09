# Ballistic Launch Angle Solver

A native solver that computes launch angles to intercept a moving target  
under gravity and quadratic drag.

Unlike vacuum or closed-form solutions, this solver works in fully nonlinear
conditions and converges even when trajectories are strongly curved.

---

## Demo (Unity)

Highly curved trajectories under strong air drag, still converging to a hit
against moving targets.

![Assets_Recordings_ballistic_demo (1)](https://github.com/user-attachments/assets/d88c2f91-09c6-48be-a950-461bfef18ed4)

Full video is available in the repository or GitHub Releases.

---

## What is this?

**ballistic-solver** computes launch angles that cause a projectile to intercept
a moving target.

The solver does not assume an analytic solution. Instead, it numerically
simulates projectile motion and iteratively adjusts the launch direction
until the closest approach to the target is minimized.

It is intended for simulations, games, and systems where robustness matters
more than closed-form elegance.

---

## How it works (high level)

1. Simulate projectile motion using RK4 integration with quadratic drag  
2. Track the closest approach between projectile and target  
3. Express the miss at closest approach as an angular residual  
4. Solve the nonlinear system using a damped least-squares method  
   (Levenberg–Marquardt with Broyden updates)  
5. Return the best-found solution, even if perfect convergence is not achieved  

Failure cases are explicitly detected and reported.

---

## Key properties

- Moving targets supported  
- Strong air resistance supported  
- No analytic assumptions  
- Robust convergence behavior  
- Explicit success / failure reporting  
- Header-only C++ core  
- Stable C ABI for multi-language use  

---

## Supported usage

The core solver is written in C++, but exposed through a stable C ABI.
This allows usage from many environments:

- C / C++
- Python (ctypes)
- C# / .NET / Unity (P/Invoke)
- Rust, Go, Node.js, Java (via FFI)

Prebuilt binaries are provided via GitHub Releases.

---

## Using prebuilt binaries

Download the archive for your platform from **Releases**.

Each release contains:

- Shared library  
  - Windows: `ballistic_c.dll`  
  - Linux: `libballistic_c.so`  
  - macOS: `libballistic_c.dylib`
- C ABI header: `ballistic_c_api.h`

The C ABI exposes a single entry point:

```c
int32_t ballistic_solve(const BallisticInputs* in, BallisticOutputs* out);
````

---

## Python usage

A minimal Python wrapper using `ctypes` is provided.

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

The native library is loaded automatically from the package.

---

## C# / Unity usage

A C# P/Invoke example is available in:

```
examples/csharp/
```

On Windows, place `ballistic_c.dll` next to the executable
(or ensure it is discoverable via PATH),
then call `ballistic_solve` via `DllImport`.

This works directly inside Unity.

---

## Build from source

```bash
cmake -S . -B build
cmake --build build -j
ctest --test-dir build
```

The shared library target is `ballistic_c`.

---

## ABI notes

* Plain C layout across the ABI boundary
* Fixed-size arrays only
* No dynamic allocation across the boundary
* The solver always returns the best-found result, even on partial failure

---

## License

MIT License
