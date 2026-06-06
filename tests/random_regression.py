import math
import random

import ballistic_solver as bs


def _dir(theta, phi):
    ct = math.cos(theta)
    return (ct * math.cos(phi), ct * math.sin(phi), math.sin(theta))


def _vacuum_pos(speed, theta, phi, t, g):
    d = _dir(theta, phi)
    return (
        speed * d[0] * t,
        speed * d[1] * t,
        speed * d[2] * t - 0.5 * g * t * t,
    )


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


def test_vacuum_stationary_matches_analytic_solution():
    params = bs.params_preset("precise")
    params.wind = (0.0, 0.0, 0.0)
    params.dt = 0.002
    params.tolMiss = 1e-3

    rel_pos = (120.0, 30.0, 5.0)
    speed = 90.0
    analytic = bs.vacuum_arc_angles_to_point(rel_pos, speed, "low", params.g)
    assert analytic["reachable"]

    result = bs.solve(rel_pos, (0.0, 0.0, 0.0), speed, 0.0, params=params)
    assert result["success"]
    assert result["miss"] <= params.tolMiss
    assert abs(result["theta"] - analytic["theta"]) < 2e-3
    assert abs(result["phi"] - analytic["phi"]) < 2e-3


def test_vacuum_moving_target_constructed_hit():
    params = bs.params_preset("precise")
    params.wind = (0.0, 0.0, 0.0)
    params.dt = 0.002
    params.tolMiss = 1e-3

    speed = 100.0
    theta = math.radians(24.0)
    phi = math.radians(18.0)
    t_hit = 1.8
    rel_vel = (3.0, -2.0, 0.5)
    hit = _vacuum_pos(speed, theta, phi, t_hit, params.g)
    rel_pos0 = tuple(hit[i] - rel_vel[i] * t_hit for i in range(3))

    result = bs.solve(rel_pos0, rel_vel, speed, 0.0, params=params)
    assert result["success"]
    assert result["miss"] <= params.tolMiss


def test_unreachable_vacuum_target_and_failure_status_examples():
    unreachable = bs.vacuum_arc_angles_to_point((10000.0, 0.0, 0.0), 50.0, "low", 9.80665)
    assert not unreachable["reachable"]

    invalid = bs.solve((100.0, 0.0, 0.0), (0.0, 0.0, 0.0), 0.0, 0.0)
    assert not invalid["success"]
    assert invalid["status"] == 1

    params = bs.params_preset("fast")
    params.tMax = 2.0
    params.maxIter = 4
    impossible = bs.solve((10000.0, 0.0, 0.0), (0.0, 0.0, 0.0), 50.0, 0.0, params=params)
    assert not impossible["success"]
    assert impossible["status"] in {2, 3, 7, 8}


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


def test_predictor_seam_matches_linear_and_accel():
    p0 = (140.0, 40.0, 8.0)
    v = (6.0, -4.0, 1.0)
    a = (0.3, 0.2, -0.1)
    v0 = 95.0
    k_drag = 0.0015

    # A constant-velocity predictor must match solve() (which assumes relVel).
    pred_cv = lambda t: (p0[0] + v[0] * t, p0[1] + v[1] * t, p0[2] + v[2] * t)
    rp = bs.solve_predicted(pred_cv, v0, k_drag)
    rs = bs.solve(p0, v, v0, k_drag)
    assert rp["success"] and rs["success"]
    assert abs(rp["theta"] - rs["theta"]) < 1e-5
    assert abs(rp["phi"] - rs["phi"]) < 1e-5

    # A constant-acceleration predictor must match solve_accel().
    pred_ca = lambda t: (
        p0[0] + v[0] * t + 0.5 * a[0] * t * t,
        p0[1] + v[1] * t + 0.5 * a[1] * t * t,
        p0[2] + v[2] * t + 0.5 * a[2] * t * t,
    )
    rp2 = bs.solve_predicted(pred_ca, v0, k_drag)
    ra = bs.solve_accel(p0, v, a, v0, k_drag)
    assert rp2["success"] and ra["success"]
    assert abs(rp2["theta"] - ra["theta"]) < 1e-6
    assert abs(rp2["phi"] - ra["phi"]) < 1e-6

    # A genuinely curved predictor must still converge and hit.
    pred_curved = lambda t: (p0[0] + v[0] * t, p0[1] + v[1] * t + 12.0 * math.sin(0.6 * t), p0[2] + v[2] * t)
    rp3 = bs.solve_predicted(pred_curved, v0, k_drag)
    assert rp3["success"]
    assert rp3["miss"] <= 1e-2


if __name__ == "__main__":
    test_vacuum_stationary_matches_analytic_solution()
    test_vacuum_moving_target_constructed_hit()
    test_unreachable_vacuum_target_and_failure_status_examples()
    test_random_linear_cases()
    test_acceleration_api_smoke()
    test_predictor_seam_matches_linear_and_accel()
    print("random_regression.py: ok")
