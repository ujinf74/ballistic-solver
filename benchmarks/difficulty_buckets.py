"""Difficulty-bucketed benchmark for the default solver (`bs.solve`).

Unlike a single aggregate number, this sweeps the genuinely hard knobs -- arc
(low / high) and drag strength (none -> heavy) -- so the easy *and* the hard
corners are both visible. Each case is constructed from a known reachable launch
(simulate a real shot, then place a moving target on it), so a sub-100% success
rate reflects the solver, not an impossible target.

Methodology (stated so the numbers can be reproduced and audited):
  * solver:     bs.solve, "precise" preset, dt=0.004, tolMiss=1e-2
  * tMax:       8.5 s for low arc, 20.0 s for high arc (plunging shots fly longer)
  * success:    result["success"] (miss <= tolMiss)
  * target:     constant-velocity moving target, reachable by construction
  * timing:     per-call time.perf_counter on the native extension, after warm-up
  * seed:       fixed (so the case set is deterministic; --dump writes it to JSON)

    python benchmarks/difficulty_buckets.py
    python benchmarks/difficulty_buckets.py --dump benchmark_cases.json

Runtime numbers depend on CPU, OS, compiler, and build type -- see the README
hardware note. The success rates and miss distributions do not.
"""

import argparse
import json
import math
import random
import statistics
import time

import ballistic_solver as bs

SEED = 20260601
N_PER_BUCKET = 1500
TOL_MISS = 1e-2

# (label, arc, k_drag); arc selects the elevation band, time-of-flight window, and
# solver branch. High-arc cases are forced past apex so the intercept is plunging.
BUCKETS = [
    ("Low arc,  no drag",       "low",  0.0),
    ("Low arc,  moderate drag", "low",  1.0e-3),
    ("Low arc,  heavy drag",    "low",  2.0e-3),
    ("High arc, no drag",       "high", 0.0),
    ("High arc, moderate drag", "high", 1.0e-3),
    ("High arc, heavy drag",    "high", 2.0e-3),
]

# elevation band and time-of-flight window per arc (matches the project's
# internal compare harness so external runs reproduce the same regimes).
ARC_THETA_DEG = {"low": (2.0, 38.0), "high": (50.0, 78.0)}
ARC_T_WINDOW = {"low": (0.8, 8.5), "high": (2.0, 12.0)}


def direction(theta, phi):
    ct = math.cos(theta)
    return (ct * math.cos(phi), ct * math.sin(phi), math.sin(theta))


def make_params(arc):
    p = bs.params_preset("precise")
    p.dt = 0.004
    p.tMax = 20.0 if arc == "high" else 8.5
    p.tolMiss = TOL_MISS
    p.arcMode = bs.ArcMode.High if arc == "high" else bs.ArcMode.Low
    return p


def make_cases(arc, k_drag, count, params, seed):
    """Construct reachable moving-target cases by simulating a real shot first.

    For the high arc the intercept time is drawn AFTER apex (apex+0.5 .. apex+8),
    so the projectile is descending onto the target -- a genuine plunging shot,
    not a steep direct-fire shot that the low branch would also solve.
    """
    rng = random.Random(seed)
    theta_min, theta_max = (math.radians(d) for d in ARC_THETA_DEG[arc])
    t_min, t_max = ARC_T_WINDOW[arc]
    cases = []
    attempts = 0
    max_attempts = max(count * 200, 1000)
    while len(cases) < count:
        attempts += 1
        if attempts > max_attempts:
            raise RuntimeError(f"could not generate {count} {arc}-arc cases within tMax={params.tMax:g}")
        speed = rng.uniform(120.0, 420.0)
        theta = rng.uniform(theta_min, theta_max)
        phi = rng.uniform(-math.pi, math.pi)
        if arc == "high":
            apex = speed * math.sin(theta) / max(params.g, 1e-9)
            t_lo = max(t_min, apex + 0.5)
            t_hi = min(min(t_max, params.tMax), apex + 8.0)
            if t_lo > t_hi:
                continue
            t_hit = rng.uniform(t_lo, t_hi)
        else:
            t_hit = rng.uniform(t_min, t_max)
        steps = max(1, round(t_hit / params.dt))
        t_hit = steps * params.dt
        rel_vel = (rng.uniform(-25.0, 25.0), rng.uniform(-25.0, 25.0), rng.uniform(-8.0, 8.0))
        d = direction(theta, phi)
        traj = bs.simulate((0.0, 0.0, 0.0), tuple(speed * x for x in d), k_drag, steps, params)
        hit = traj["r"][-1]
        rel_pos0 = tuple(hit[i] - rel_vel[i] * t_hit for i in range(3))
        cases.append((rel_pos0, rel_vel, speed, k_drag))
    return cases


def percentile(values, q):
    values = sorted(values)
    return values[min(len(values) - 1, round((len(values) - 1) * q))]


def run_bucket(label, arc, k_drag):
    params = make_params(arc)
    cases = make_cases(arc, k_drag, N_PER_BUCKET, params, SEED)

    for rel_pos0, rel_vel, speed, kd in cases[:5]:
        bs.solve(rel_pos0, rel_vel, speed, kd, arcMode=arc, params=params)

    times, misses, successes = [], [], 0
    for rel_pos0, rel_vel, speed, kd in cases:
        t0 = time.perf_counter()
        result = bs.solve(rel_pos0, rel_vel, speed, kd, arcMode=arc, params=params)
        times.append(time.perf_counter() - t0)
        misses.append(result["miss"])
        successes += bool(result["success"])

    return {
        "label": label,
        "n": len(cases),
        "success": successes,
        "median_ms": statistics.median(times) * 1000.0,
        "p95_ms": percentile(times, 0.95) * 1000.0,
        "p99_ms": percentile(times, 0.99) * 1000.0,
        "p95_miss": percentile(misses, 0.95),
        "cases": cases,
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dump", metavar="PATH", help="write the fixed case set to a JSON file")
    args = ap.parse_args()

    rows = [run_bucket(*b) for b in BUCKETS]

    print(f"Default solver (bs.solve), precise preset, dt=0.004, tolMiss={TOL_MISS:g}, seed={SEED}")
    print(f"tMax = 8.5 s (low arc) / 20.0 s (high arc); {N_PER_BUCKET} reachable moving-target cases per bucket.\n")
    print("| Bucket | Success | Median | P95 | P99 | P95 miss |")
    print("|---|---:|---:|---:|---:|---:|")
    for r in rows:
        pct = 100.0 * r["success"] / r["n"]
        print(
            f"| {r['label']} | {r['success']}/{r['n']} ({pct:.2f}%) "
            f"| {r['median_ms']:.3f} ms | {r['p95_ms']:.3f} ms | {r['p99_ms']:.3f} ms "
            f"| {r['p95_miss']:.3e} m |"
        )

    if args.dump:
        payload = {
            "seed": SEED,
            "tolMiss": TOL_MISS,
            "params": {"preset": "precise", "dt": 0.004, "tMax": 8.0},
            "buckets": [
                {
                    "label": r["label"],
                    "cases": [
                        {"relPos0": rp, "relVel": rv, "v0": v0, "kDrag": kd}
                        for rp, rv, v0, kd in r["cases"]
                    ],
                }
                for r in rows
            ],
        }
        with open(args.dump, "w", encoding="utf-8") as f:
            json.dump(payload, f)
        total = sum(r["n"] for r in rows)
        print(f"\nWrote {total} cases to {args.dump}")


if __name__ == "__main__":
    main()
