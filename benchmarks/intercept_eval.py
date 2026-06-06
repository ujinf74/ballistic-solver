"""CIWS / gun-AA engagement realism check.

Kinematic feasibility test (NOT a validated weapon model): can a position-only
radar track + the predictor seam produce a firing solution that actually
intercepts a fast, possibly maneuvering target?

Pipeline per scenario:
  1. Ground-truth target trajectory (gun at origin).
  2. Noisy position-only radar measurements at 100 Hz -> bs.TargetTracker.
  3. At fire time (target inside engagement range) compute three firing
     solutions: ideal (true pos/vel/acc), tracker (from noisy history),
     and a naive constant-velocity baseline (two-point finite difference).
  4. Simulate the actual projectile (RK4) and report its closest approach to
     the TRUE (continuing) target -- the real miss distance.

    python benchmarks/ciws_eval.py

Assumptions (stated, not authoritative): 20 mm-class round, muzzle 1050 m/s,
quadratic drag from a 0.30 Cd / 0.10 kg / 20 mm projectile, sea-level air.
"""

import math
import random
import time

import ballistic_solver as bs


V0 = 1050.0
K_DRAG = bs.k_drag_from_physical(airDensity=1.225, dragCoefficient=0.30, area=math.pi * 0.010 ** 2, mass=0.10)
G = 9.80665


def direction(theta, phi):
    ct = math.cos(theta)
    return (ct * math.cos(phi), ct * math.sin(phi), math.sin(theta))


def solve_params():
    p = bs.params_preset("precise")
    p.dt = 0.004
    p.tMax = 6.0
    return p


def real_miss(theta, phi, true_at, dt=0.002, tmax=6.0):
    d = direction(theta, phi)
    vv = (V0 * d[0], V0 * d[1], V0 * d[2])
    p = bs.params_preset("balanced")
    p.dt = dt
    p.g = G
    steps = int(tmax / dt)
    traj = bs.simulate((0.0, 0.0, 0.0), vv, K_DRAG, steps, p)
    R, T = traj["r"], traj["t"]
    best = float("inf")
    for i in range(len(R)):
        best = min(best, math.dist(R[i], true_at(T[i])))
    return best


# ----- scenarios: pos/vel/acc as functions of global time (gun-frame) -----
def straight_inbound():
    s, v = (2600.0, 200.0, 12.0), (-296.0, -23.0, 0.0)
    return {
        "name": "sea-skimmer, straight",
        "pos": lambda t: (s[0] + v[0] * t, s[1] + v[1] * t, s[2] + v[2] * t),
        "vel": lambda t: v,
        "acc": lambda t: (0.0, 0.0, 0.0),
        "fire_range": 1600.0,
        "tmax": 6.0,
    }


def weaving_inbound(amp=25.0, omega=2.5):
    s, v = (2600.0, 200.0, 12.0), (-296.0, -23.0, 0.0)
    return {
        "name": "sea-skimmer, weaving ~16g",
        "pos": lambda t: (s[0] + v[0] * t, s[1] + v[1] * t + amp * math.sin(omega * t), s[2]),
        "vel": lambda t: (v[0], v[1] + amp * omega * math.cos(omega * t), 0.0),
        "acc": lambda t: (0.0, -amp * omega * omega * math.sin(omega * t), 0.0),
        "fire_range": 1600.0,
        "tmax": 6.0,
    }


def crossing():
    s, v = (1500.0, -1500.0, 300.0), (-150.0, 250.0, 0.0)
    return {
        "name": "aircraft, crossing",
        "pos": lambda t: (s[0] + v[0] * t, s[1] + v[1] * t, s[2] + v[2] * t),
        "vel": lambda t: v,
        "acc": lambda t: (0.0, 0.0, 0.0),
        "fire_range": 1500.0,
        "tmax": 8.0,
    }


def engage(scn, noise, q, seed):
    rng = random.Random(seed)
    tracker = bs.TargetTracker(processNoise=q, measNoise=max(noise * noise, 1e-6))
    P = solve_params()

    dt_meas = 0.01
    prev = prev_t = None
    steps = int(scn["tmax"] / dt_meas)
    for k in range(steps):
        t = k * dt_meas
        true = scn["pos"](t)
        meas = tuple(true[i] + rng.gauss(0.0, noise) for i in range(3))
        tracker.update(t, meas)

        if math.dist(true, (0.0, 0.0, 0.0)) < scn["fire_range"]:
            ideal = bs.solve_accel(scn["pos"](t), scn["vel"](t), scn["acc"](t), V0, K_DRAG, params=P)
            t0 = time.perf_counter()
            trk = tracker.solve(V0, K_DRAG, params=P)
            solve_ms = (time.perf_counter() - t0) * 1000.0
            cv = None
            if prev is not None:
                dtm = t - prev_t
                vcv = tuple((meas[i] - prev[i]) / dtm for i in range(3))
                cv = bs.solve(meas, vcv, V0, K_DRAG, params=P)

            true_at = lambda s: scn["pos"](t + s)
            return {
                "tStar": trk["tStar"],
                "ideal": real_miss(ideal["theta"], ideal["phi"], true_at),
                "tracker": real_miss(trk["theta"], trk["phi"], true_at),
                "cv": real_miss(cv["theta"], cv["phi"], true_at) if cv else float("nan"),
                "solve_ms": solve_ms,
            }
        prev, prev_t = meas, t
    return None


def main():
    print(f"v0={V0:.0f} m/s, kDrag={K_DRAG:.3e}; real miss = projectile closest approach to TRUE target")
    print(f"{'scenario':<26} {'noise':>6} {'tof s':>6} {'ideal m':>9} {'tracker m':>10} {'cv m':>9} {'solve ms':>9}")
    for scn_fn in (straight_inbound, weaving_inbound, crossing):
        for noise in (1.0, 3.0):
            scn = scn_fn()
            # NOTE: process noise q must be matched to target dynamics. Low q
            # (~1) denoises non-maneuvering tracks to ~1-2 m aim error; higher q
            # tracks maneuvers but degrades steady accuracy. No fixed q is
            # optimal for all -- see the q-sweep discussion / motivates IMM.
            r = engage(scn, noise=noise, q=5.0, seed=11)
            if r is None:
                print(f"{scn['name']:<26} {noise:>6.1f}  (never entered fire range)")
                continue
            print(
                f"{scn['name']:<26} {noise:>6.1f} {r['tStar']:>6.2f} "
                f"{r['ideal']:>9.2f} {r['tracker']:>10.2f} {r['cv']:>9.2f} {r['solve_ms']:>9.3f}"
            )
    print("note: 'ideal' uses the true instantaneous pos/vel/acc; a large 'ideal' miss")
    print("means the constant-acceleration prediction cannot follow the maneuver over")
    print("the time of flight -- a fundamental limit of second-order lead, not a solver bug.")


if __name__ == "__main__":
    main()
