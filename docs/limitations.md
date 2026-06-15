# Limitations

`ballistic-solver` is a point-mass intercept solver for games, simulation,
robotics, and research. It is deliberately scoped, and the assumptions below are
part of that scope rather than bugs. Knowing them up front is the best way to
decide whether it fits your use case.

## Physical model

- **Point mass only.** The projectile is a point with no orientation, no rigid-
  body dynamics, and no lift. There is **no spin drift / Magnus effect**.
- **Quadratic drag with a single scalar coefficient** `kDrag`. There is **no
  Mach-dependent drag table** and no per-segment drag curve; transonic/supersonic
  rounds whose `Cd` varies strongly with speed are not modelled accurately.
- **Constant gravity.** No altitude variation, no Coriolis term.
- **Uniform, constant wind.** Wind is a single air-velocity vector; there is no
  gust, shear, or altitude profile.
- **Air density is folded into `kDrag`** and is constant for the flight; it does
  not vary with altitude during the trajectory.

## Solve scope

- **Two unknowns: launch angles `(θ, φ)` only.** Muzzle speed `v0` is an input,
  not solved for; the solver does not pick a charge/speed to make a target
  reachable.
- **Target motion is constant-velocity or constant-acceleration** (or whatever a
  user-supplied predictor returns via `solve_predicted`). A freely maneuvering
  target is only handled to the extent the predictor models it.
- **One arc branch at a time.** With `arcMode` fixed, the solver finds a solution
  on that branch; a solution that exists only on the other branch is not returned
  unless you switch `arcMode`.
- **No terrain or obstacle collision.** Trajectories are not clipped against
  ground or geometry; an angle that "hits" may pass through terrain in your
  world.

## Numerical behaviour

- **Fixed-step RK4.** Accuracy depends on `dt`, `tMax`, and `tolMiss`. Too large
  a `dt` or too short a `tMax` can miss or fail to converge.
- **Best-effort, not guaranteed.** For unreachable or near-limit targets the
  solver returns its closest iterate with `success = false`; it does not prove
  that no solution exists.
- **Your runtime must match the solver's physics.** If the consuming engine
  integrates with a different scheme or timestep, the computed angles may not hit
  even when the solver reports success. See the "Solver validity note" in the
  README.

## Not a fire-control system

This is a numerical library for simulated and robotic targeting and for research.
It is **not** a certified or safety-rated fire-control system and is not intended
for use as one.
