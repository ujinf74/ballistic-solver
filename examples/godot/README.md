# Godot 4 GDExtension Example

This example exposes the C ABI solver to Godot 4 through a small GDExtension
wrapper. It is intentionally kept outside the main CMake build so the core
library does not depend on Godot or `godot-cpp`.

The wrapper accepts and returns Godot coordinates: `x` right/forward range,
`y` up, and `z` lateral. Internally it converts to the core solver's
`x/y` horizontal plane and `z` up convention.

## Layout

```text
examples/godot/
  project.godot
  addons/ballistic_solver/ballistic_solver.gdextension
  addons/ballistic_solver/plugin.cfg
  demo/ballistic_demo.gd     # perfect-info lead-fire demo
  demo/intercept_demo.gd     # noisy-track tracker intercept demo (default scene)
  src/ballistic_solver_gd.h
  src/ballistic_solver_gd.cpp
  SConstruct
```

## Classes

* `BallisticSolver` — `solve(...)` and `simulate_from_angles(...)` over the C ABI.
* `BallisticTracker` — a position-only constant-acceleration tracker wrapping the
  C++ core directly (`TargetTracker` from `bs/kalman.hpp`): `configure(process,
  meas)`, `update(t, rel_pos)`, `predict(tau)`, `estimated_position/velocity/
  acceleration()`, and `solve(v0, k_drag, arc_mode)` for lead-fire via the
  predictor seam. Adding it changed the extension source, so **rebuild the
  GDExtension** (see below) before the intercept demo will run.

The intercept demo (`demo/intercept_demo.tscn`, the default main scene) feeds the
tracker NOISY position-only measurements of a moving target and fires a
continuous burst from the tracker's lead solution; rounds are scored against the
true target, so green markers reflect real intercepts and the HUD shows the
running hit rate. The same machinery applies across domains (robotics, perimeter
security, agricultural deterrence, games/sim, defense).

## Build outline

1. Build or download the native `ballistic_solver` shared library for your
   platform.
2. Build `godot-cpp` for the same Godot version and target.
3. From this directory, set `GODOT_CPP` to the `godot-cpp` checkout and run:

```bash
scons platform=windows target=template_debug
```

PowerShell example:

```powershell
$env:GODOT_CPP = "C:\path\to\godot-cpp"
$env:BALLISTIC_SOLVER_BUILD_DIR = "..\..\build\Release"
python -m SCons platform=windows target=template_debug arch=x86_64
```

Use `platform=linux` or `platform=macos` on other systems. If the native solver
was built outside `../../build/Release`, set `BALLISTIC_SOLVER_BUILD_DIR` to the
directory containing `ballistic_solver.lib` before running `scons`.

Copy the `ballistic_solver` shared library next to the generated GDExtension
library. The `.gdextension` manifest declares this runtime dependency so Godot
loads the native solver before loading the wrapper:

```bash
copy ..\..\build\Release\ballistic_solver.dll addons\ballistic_solver\bin\
```

Then open this folder as a Godot project.

## Smoke test

After building, the first run must import the project so Godot registers the
GDExtension; then run the main scene headless (`--quit-after` counts frames):

```bash
godot --headless --editor --path . --quit       # one-time import (registers the extension)
godot --headless --path . --quit-after 250       # run the main scene
```

Verified with Godot 4.2.2 (headless): the extension loads, `BallisticTracker`
and `BallisticSolver` resolve, and the scene runs without script/extension
errors (headless emits harmless dummy-renderer "Parameter m is null" mesh
warnings). For the visuals, open the folder in the Godot editor and press Play.

The intercept scene tracks a moving target from noisy measurements and fires a burst;
expect red radar dots around the target, a yellow predicted-intercept marker,
grey rounds in flight, and green hit markers with a hit-rate HUD. The simpler
`res://demo/ballistic_demo.tscn` instantiates `BallisticSolver`, calls `solve`,
draws the sampled drag trajectory, and prints the returned dictionary; a valid
run shows `"success": true` and a cyan trajectory, yellow lead line, red target,
green intercept marker, and blue projectile.

The wrapper returns a `Dictionary` with `success`, `status`, `theta`, `phi`,
`miss`, `t_star`, `iterations`, `accepted_steps`, `last_lambda`, `last_alpha`,
and `message`.

## Addon use

For another Godot project, copy `addons/ballistic_solver/` after building the
native libraries. The demo scene is separate and is not required at runtime.
