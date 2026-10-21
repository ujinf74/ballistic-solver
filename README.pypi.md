# ballistic-solver

[![PyPI](https://img.shields.io/pypi/v/ballistic-solver)](https://pypi.org/project/ballistic-solver/)
[![CI](https://github.com/ujinf74/ballistic-solver/actions/workflows/ci.yml/badge.svg?branch=main)](https://github.com/ujinf74/ballistic-solver/actions/workflows/ci.yml)

Given a moving target, a muzzle speed and a drag coefficient, this returns the
launch angles that hit it. It integrates the projectile and solves the intercept
numerically, so it stays accurate where closed-form vacuum formulas drift: strong
drag, high arcs, wind, and targets that keep moving while the round is in flight.

The core is C++17 with a stable C ABI; this package is its native Python
extension, not a pure-Python reimplementation.

## Install

```bash
pip install ballistic-solver
```

Wheels cover CPython 3.10-3.15 on Linux (`x86_64`, `aarch64`; manylinux and
musllinux), macOS (`universal2`) and Windows (`AMD64`, `ARM64`). Other
combinations build from the sdist, which needs CMake and a C++17 compiler.

## Use

```python
import ballistic_solver as bs

result = bs.solve(
    relPos0=(120, 30, 5),    # target position relative to the muzzle, metres
    relVel=(2, -1, 0),       # target velocity, m/s
    v0=90,                   # muzzle speed, m/s
    kDrag=0.002,             # quadratic drag coefficient
)

print(result["success"], result["theta"], result["phi"])
```

`theta` and `phi` are elevation and azimuth in radians. `miss` is the
closest-approach distance and `tStar` the time it occurs; `success` reports
whether the miss landed inside the tolerance. A best-effort answer comes back
even when it does not, so an unreachable target degrades rather than throws.

Derive the drag coefficient from physical quantities if you have them:

```python
kDrag = bs.k_drag_from_physical(
    airDensity=1.225, dragCoefficient=0.30, area=0.00426, mass=0.145,
)
```

Three presets trade speed against tolerance — `bs.params_preset("fast")`,
`"balanced"`, `"precise"` — and `bs.solve_accel` accepts a constant target
acceleration. For a noisy position-only feed, `bs.TargetTracker` denoises the
track and fires through the same solver without you supplying velocity.

## Scope

A point-mass model: quadratic drag with a single coefficient, constant gravity,
uniform wind. No spin drift or Magnus effect, no Mach-dependent drag table, no
terrain collision. It solves launch angles only — muzzle speed and arc branch
are inputs. Intended for games, simulation, robotics and research; it is not a
certified fire-control system.

## Documentation

The full API reference, the C ABI, the C++ surface, the Unity and Godot
bindings, benchmarks and the numerical method all live in the repository:

**https://github.com/ujinf74/ballistic-solver**

- [Numerical method](https://github.com/ujinf74/ballistic-solver/blob/main/docs/numerical_method.md)
- [Benchmarks](https://github.com/ujinf74/ballistic-solver/blob/main/docs/benchmarks.md)
- [Limitations](https://github.com/ujinf74/ballistic-solver/blob/main/docs/limitations.md)
- [Changelog](https://github.com/ujinf74/ballistic-solver/blob/main/CHANGELOG.md)

## License

MIT
