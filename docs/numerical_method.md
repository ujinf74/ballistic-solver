# Numerical method

This document describes how the default solver (`bs.solve`) computes launch
angles. It is meant to be read alongside the source in
`include/bs/solve_coord_lead.hpp`.

## 1. Problem

A projectile is fired from the origin at a fixed muzzle speed `v0` in a direction
`d(θ, φ)`, where `θ` is elevation and `φ` is azimuth. It flies under constant
gravity `g`, quadratic air drag with coefficient `kDrag`, and an optional
constant wind. The target moves with a known relative motion

```
aim(t) = relPos0 + relVel · t + ½ · relAcc · t²
```

(`relAcc` is zero for the constant-velocity API). The task is to choose `(θ, φ)`
so the projectile passes as close as possible to the target.

The two unknowns are the launch angles only. Muzzle speed and time of flight are
not free variables: time of flight is whatever the simulation produces, and the
intercept instant is found as the point of closest approach (below).

## 2. Residual

For a candidate `(θ, φ)` the projectile trajectory `r_p(t)` is integrated with a
fixed-step 4th-order Runge–Kutta scheme:

```
a(v) = (0, 0, −g) − kDrag · |v − wind| · (v − wind)
```

The relative miss is `m(t) = r_p(t) − aim(t)`. The solver tracks the time of
closest approach `t* = argmin_t ||m(t)||` and defines the **residual** as the 3D
miss vector at that instant:

```
r = m(t*) ∈ ℝ³,      miss = ||r||
```

Using the coordinate miss vector directly as the residual is the key choice: it
is well-conditioned and varies smoothly with the angles, unlike an
angle-space residual that has to be reconstructed through an auxiliary inverse.

The solve is therefore a 2-input / 3-output root problem `r(θ, φ) → 0`, handled
by a quasi-Newton iteration with a 2×3 inverse-Jacobian operator `M`.

## 3. Vacuum seed

The iteration starts from a closed-form vacuum lead
(`initial_guess_vacuum_lead_acc`): the moving-target quartic gives the launch
angles that would intercept in the absence of drag. `θ` is clamped to
`[thetaMin, thetaMax]` and `φ` is wrapped to `(−π, π]`. This places the warm
start near the true solution for low-to-moderate drag, which is where most of
the speed advantage comes from.

## 4. Inverse-Jacobian preconditioning (vacuum-analytic seed)

Rather than build a Jacobian by finite-differencing the *simulated* trajectory
(which costs extra RK4 integrations per column), the initial inverse Jacobian

```
M = ∂(θ, φ) / ∂(aim)        (2×3)
```

is finite-differenced from the **closed-form** vacuum map
`vacuum_arc_angles_to_point` at the current aim point `aim(t*)`. No trajectory
integration is used to form `M`. `M` answers "if the aim point moves by δ in
world space, how do the launch angles change", and is used to turn a coordinate
miss back into an angle correction. It is an analytic surrogate / preconditioner,
not the exact Jacobian of the drag dynamics.

## 5. Quasi-Newton step and Broyden update

Each iteration takes a **full step** (no Levenberg–Marquardt damping, no line
search):

```
δθ = −(M[0] · r)
δφ = −(M[1] · r)
θ ← clamp(θ + δθ, thetaMin, thetaMax)
φ ← wrap(φ + δφ)
```

After re-evaluating the residual, `M` is refined with a **bad-Broyden inverse
update** so it tracks the true coordinate map (the vacuum surrogate is only the
seed):

```
s = Δ(θ, φ)        (2-vector)
y = Δr             (3-vector)
M ← M + ((s − M·y) · yᵀ) / (yᵀ·y)        when yᵀ·y is non-negligible
```

The best (lowest-miss) iterate seen so far is always retained.

## 6. Stopping condition

The iteration stops and reports success (`status = Ok`) as soon as
`miss ≤ tolMiss`. Otherwise it runs up to `maxIter` iterations and returns the
best-so-far iterate as a best-effort result. If the vacuum inverse Jacobian
cannot be formed at the seed, it returns `JacobianFailed`; if the iteration
limit is hit without reaching tolerance, `MaxIterReached`.

`tolMiss`, `maxIter`, `dt`, and `tMax` come from the chosen preset
(`fast` / `balanced` / `precise`) or from an explicit `BallisticParams`.

## 7. Multistart fallback

If the warm start does not converge, the solver retries from an
arc-appropriate elevation grid (azimuth taken from the vacuum seed):

- High arc: `{50, 54, 58, 62, 66, 70, 74, 78}°`
- Low arc:  `{6, 12, 18, 24, 30, 36}°`

Each grid point runs the same single solve; the first success is returned,
otherwise the lowest-miss result across all attempts. This grid is what gives
the method its robustness on hard cases (strong drag, near-limit geometry) while
keeping the common case to a single warm-started solve.

## 8. Known failure modes

- Invalid inputs (`v0 ≤ 0`, non-positive `g` / `dt` / `tMax` / `maxIter`, or an
  empty `[thetaMin, thetaMax]` range) → `InvalidInput`.
- Geometrically unreachable vacuum targets (the seed map has no solution).
- Targets whose intercept lies beyond `tMax`.
- Strong-drag or high-arc cases where the requested arc branch cannot be
  maintained.
- A tolerance that is tighter than the iteration budget can reach.

In all of these the solver still returns its best-effort iterate; callers should
check `success`, `status`, `miss`, and `message` together. See
[limitations.md](limitations.md) for modelling assumptions.

## 9. Relation to the auxiliary-residual method (`solve_aux`)

The previous default (`bs.solve_aux`) solves an **auxiliary-residual**
formulation: the coordinate miss is routed through the vacuum inverse to form an
angle-space residual, which is solved by damped least squares
(Levenberg–Marquardt) with the same multistart fallback. The coordinate-residual
core above instead uses the raw 3D miss as the residual and the vacuum map only
as seed and preconditioner. At equal robustness it is simpler and typically
~15–25% faster on the benchmark cases. `solve_aux` is kept for compatibility and
for reproducing the auxiliary-residual results.
