import math
import random
import statistics
import time

import ballistic_solver as bs


def direction(theta, phi):
    ct = math.cos(theta)
    return (ct * math.cos(phi), ct * math.sin(phi), math.sin(theta))


def make_cases(count, params, seed=20260502):
    rng = random.Random(seed)
    cases = []

    for _ in range(count):
        speed = rng.uniform(120.0, 420.0)
        theta = rng.uniform(math.radians(2.0), math.radians(38.0))
        phi = rng.uniform(-math.pi, math.pi)
        t_hit = rng.uniform(0.8, 8.5)
        steps = max(1, round(t_hit / params.dt))
        t_hit = steps * params.dt

        rel_vel = (rng.uniform(-25.0, 25.0), rng.uniform(-25.0, 25.0), rng.uniform(-8.0, 8.0))
        k_drag = rng.choice([0.0, 0.001, 0.002])
        d = direction(theta, phi)
        traj = bs.simulate((0.0, 0.0, 0.0), tuple(speed * x for x in d), k_drag, steps, params)
        hit = traj["r"][-1]
        rel_pos0 = tuple(hit[i] - rel_vel[i] * t_hit for i in range(3))
        cases.append((rel_pos0, rel_vel, speed, k_drag))

    return cases


def percentile(values, q):
    values = sorted(values)
    return values[min(len(values) - 1, round((len(values) - 1) * q))]


def run(preset, count):
    params = bs.params_preset(preset)
    cases = make_cases(count, params)
    times = []
    misses = []
    successes = 0

    for rel_pos0, rel_vel, speed, k_drag in cases[:5]:
        bs.solve(rel_pos0, rel_vel, speed, k_drag, params=params)

    for rel_pos0, rel_vel, speed, k_drag in cases:
        t0 = time.perf_counter()
        result = bs.solve(rel_pos0, rel_vel, speed, k_drag, params=params)
        times.append(time.perf_counter() - t0)
        misses.append(result["miss"])
        successes += bool(result["success"])

    print(
        f"{preset}: success={successes}/{count} "
        f"median={statistics.median(times) * 1000:.3f}ms "
        f"p95={percentile(times, 0.95) * 1000:.3f}ms "
        f"median_miss={statistics.median(misses):.3e} "
        f"p95_miss={percentile(misses, 0.95):.3e}"
    )


if __name__ == "__main__":
    for name in ("fast", "balanced", "precise"):
        run(name, 500)
