"""Moving-target intercept realism check (general purpose).

Kinematic feasibility test for lead-firing at a moving target from a noisy,
position-only track. The same solver + tracker applies across domains -- it only
sees relative kinematics -- so the scenarios below span very different speed and
range regimes:

  * robotics     -- a slow ground vehicle (e.g. capture/net launcher)
  * security     -- a crossing perimeter intruder
  * agriculture  -- a maneuvering field drone (bird/pest deterrence)
  * defense      -- a fast, hard-weaving inbound threat (CIWS-like)

Pipeline per scenario:
  1. Ground-truth target trajectory (launcher at origin).
  2. Noisy position-only measurements at 100 Hz -> bs.TargetTracker.
  3. At fire time (target inside engagement range) compute three firing
     solutions: ideal (true pos/vel/acc), tracker (from noisy history), and a
     naive constant-velocity baseline (two-point finite difference).
  4. Simulate the actual projectile (RK4) and report its closest approach to the
     TRUE (continuing) target -- the real miss distance.

    python benchmarks/intercept_eval.py

These are kinematic feasibility numbers, not validated domain models.
"""

import math
import random
import time

import ballistic_solver as bs

G = 9.80665
ROUND_20MM = bs.k_drag_from_physical(airDensity=1.225, dragCoefficient=0.30, area=math.pi * 0.010 ** 2, mass=0.10)


def direction(theta, phi):
    ct = math.cos(theta)
    return (ct * math.cos(phi), ct * math.sin(phi), math.sin(theta))


def solve_params():
    p = bs.params_preset("precise")
    p.dt = 0.004
    p.tMax = 8.0
    return p


def real_miss(v0, k_drag, theta, phi, true_at, dt=0.002, tmax=8.0):
    d = direction(theta, phi)
    vv = (v0 * d[0], v0 * d[1], v0 * d[2])
    p = bs.params_preset("balanced")
    p.dt = dt
    p.g = G
    steps = int(tmax / dt)
    traj = bs.simulate((0.0, 0.0, 0.0), vv, k_drag, steps, p)
    R, T = traj["r"], traj["t"]
    best = float("inf")
    for i in range(len(R)):
        best = min(best, math.dist(R[i], true_at(T[i])))
    return best


# ----- scenarios (positions in solver frame: x,y horizontal, z up) -----
def ground_robot():
    s, v, a = (90.0, 30.0, 1.5), (-8.0, -3.0, 0.0), (0.0, 0.4, 0.0)
    return {
        "name": "ground vehicle", "domain": "robotics", "v0": 60.0, "k_drag": 0.010,
        "pos": lambda t: (s[0] + v[0] * t, s[1] + v[1] * t + 0.5 * a[1] * t * t, s[2]),
        "vel": lambda t: (v[0], v[1] + a[1] * t, 0.0),
        "acc": lambda t: a,
        "fire_range": 85.0, "tmax": 8.0,
    }


def perimeter_crossing():
    s, v = (60.0, -150.0, 2.0), (-6.0, 13.0, 0.0)
    return {
        "name": "crossing intruder", "domain": "security", "v0": 120.0, "k_drag": 0.006,
        "pos": lambda t: (s[0] + v[0] * t, s[1] + v[1] * t, s[2]),
        "vel": lambda t: v,
        "acc": lambda t: (0.0, 0.0, 0.0),
        "fire_range": 150.0, "tmax": 8.0,
    }


def field_drone(amp=6.0, omega=1.5):
    s, v = (150.0, 40.0, 25.0), (-16.0, -8.0, 0.0)
    return {
        "name": "maneuvering drone", "domain": "agriculture", "v0": 110.0, "k_drag": 0.008,
        "pos": lambda t: (s[0] + v[0] * t, s[1] + v[1] * t + amp * math.sin(omega * t), s[2]),
        "vel": lambda t: (v[0], v[1] + amp * omega * math.cos(omega * t), 0.0),
        "acc": lambda t: (0.0, -amp * omega * omega * math.sin(omega * t), 0.0),
        "fire_range": 150.0, "tmax": 8.0,
    }


def inbound_weaving(amp=25.0, omega=2.5):
    s, v = (2600.0, 200.0, 12.0), (-296.0, -23.0, 0.0)
    return {
        "name": "fast inbound, weaving", "domain": "defense", "v0": 1050.0, "k_drag": ROUND_20MM,
        "pos": lambda t: (s[0] + v[0] * t, s[1] + v[1] * t + amp * math.sin(omega * t), s[2]),
        "vel": lambda t: (v[0], v[1] + amp * omega * math.cos(omega * t), 0.0),
        "acc": lambda t: (0.0, -amp * omega * omega * math.sin(omega * t), 0.0),
        "fire_range": 1600.0, "tmax": 6.0,
    }


SCENARIOS = [ground_robot, perimeter_crossing, field_drone, inbound_weaving]


def engage(scn, noise, q, seed):
    rng = random.Random(seed)
    tracker = bs.TargetTracker(processNoise=q, measNoise=max(noise * noise, 1e-6))
    P = solve_params()
    v0, k_drag = scn["v0"], scn["k_drag"]

    dt_meas = 0.01
    prev = prev_t = None
    for k in range(int(scn["tmax"] / dt_meas)):
        t = k * dt_meas
        true = scn["pos"](t)
        meas = tuple(true[i] + rng.gauss(0.0, noise) for i in range(3))
        tracker.update(t, meas)

        if math.dist(true, (0.0, 0.0, 0.0)) < scn["fire_range"]:
            true_at = lambda s, t=t: scn["pos"](t + s)
            ideal = bs.solve_accel(scn["pos"](t), scn["vel"](t), scn["acc"](t), v0, k_drag, params=P)
            t0 = time.perf_counter()
            trk = tracker.solve(v0, k_drag, params=P)
            solve_ms = (time.perf_counter() - t0) * 1000.0
            cv = None
            if prev is not None:
                dtm = t - prev_t
                vcv = tuple((meas[i] - prev[i]) / dtm for i in range(3))
                cv = bs.solve(meas, vcv, v0, k_drag, params=P)
            return {
                "tStar": trk["tStar"],
                "ideal": real_miss(v0, k_drag, ideal["theta"], ideal["phi"], true_at),
                "tracker": real_miss(v0, k_drag, trk["theta"], trk["phi"], true_at),
                "cv": real_miss(v0, k_drag, cv["theta"], cv["phi"], true_at) if cv else float("nan"),
                "solve_ms": solve_ms,
            }
        prev, prev_t = meas, t
    return None


def main():
    print("real miss = projectile closest approach to the TRUE target (lower is better)")
    print(f"{'scenario':<22} {'domain':<12} {'v0':>6} {'tof s':>6} {'ideal m':>9} {'tracker m':>10} {'cv m':>9} {'ms':>6}")
    for scn_fn in SCENARIOS:
        for noise in (1.0, 3.0):
            scn = scn_fn()
            # process noise q is tuned to the dynamics; low for steady, high for maneuver.
            r = engage(scn, noise=noise, q=5.0, seed=11)
            if r is None:
                print(f"{scn['name']:<22} {scn['domain']:<12}  (never entered fire range)")
                continue
            print(
                f"{scn['name']:<22} {scn['domain']:<12} {scn['v0']:>6.0f} {r['tStar']:>6.2f} "
                f"{r['ideal']:>9.2f} {r['tracker']:>10.2f} {r['cv']:>9.2f} {r['solve_ms']:>6.2f}  noise={noise:.0f}m"
            )
    print("note: a large 'ideal' miss (true state) means the constant-acceleration prediction")
    print("cannot follow the maneuver over the time of flight -- a fundamental limit of")
    print("second-order lead, not a solver bug. q must be matched to the target dynamics.")


if __name__ == "__main__":
    main()
