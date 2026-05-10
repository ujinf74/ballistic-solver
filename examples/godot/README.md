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
  demo/ballistic_demo.gd
  src/ballistic_solver_gd.h
  src/ballistic_solver_gd.cpp
  SConstruct
```

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

After building, run the demo scene from the command line:

```bash
godot --path . res://demo/ballistic_demo.tscn --quit-after 3
```

The scene instantiates `BallisticSolver`, calls `solve`, draws the sampled drag
trajectory, and prints the returned dictionary to the Godot output log. A valid
run should include `"success": true` and a visible cyan trajectory, yellow lead
line, red target, green intercept marker, and blue projectile.

The wrapper returns a `Dictionary` with `success`, `status`, `theta`, `phi`,
`miss`, `t_star`, `iterations`, `accepted_steps`, `last_lambda`, `last_alpha`,
and `message`.

## Addon use

For another Godot project, copy `addons/ballistic_solver/` after building the
native libraries. The demo scene is separate and is not required at runtime.
