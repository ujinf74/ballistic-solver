import math
import random

import ballistic_solver as bs


def _dir(theta, phi):
    ct = math.cos(theta)
    return (ct * math.cos(phi), ct * math.sin(phi), math.sin(theta))


def _case(rng):
    speed = rng.uniform(80.0, 180.0)
    theta = rng.uniform(math.radians(3.0), math.radians(35.0))
    phi = rng.uniform(-math.pi, math.pi)
    t_hit = rng.uniform(0.8, 4.0)
    rel_vel = (rng.uniform(-8.0, 8.0), rng.uniform(-8.0, 8.0), rng.uniform(-2.0, 2.0))
    k_drag = rng.choice([0.0, 0.0005, 0.001])

    p = bs.params_preset("precise")
    p.dt = 0.005
    steps = max(1, round(t_hit / p.dt))
    t_hit = steps * p.dt

    d = _dir(theta, phi)
    traj = bs.simulate((0.0, 0.0, 0.0), tuple(speed * x for x in d), k_drag, steps, p)
    hit = traj["r"][-1]
    rel_pos0 = tuple(hit[i] - rel_vel[i] * t_hit for i in range(3))
    return rel_pos0, rel_vel, speed, k_drag, p


def test_random_linear_cases():
    rng = random.Random(20260502)
    total = 200
    successes = 0
    misses = []

    for _ in range(total):
        rel_pos0, rel_vel, speed, k_drag, params = _case(rng)
        result = bs.solve(rel_pos0, rel_vel, speed, k_drag, params=params)
        successes += bool(result["success"])
        misses.append(result["miss"])

    assert successes / total >= 0.98
    assert sorted(misses)[int(total * 0.95)] < 1e-3


def test_acceleration_api_smoke():
    params = bs.params_preset("precise")
    result = bs.solve_accel(
        relPos0=(120.0, 30.0, 5.0),
        relVel=(2.0, -1.0, 0.0),
        relAcc=(0.0, 0.2, 0.0),
        v0=90.0,
        kDrag=0.002,
        params=params,
    )

    assert math.isfinite(result["theta"])
    assert math.isfinite(result["phi"])
    assert math.isfinite(result["miss"])
