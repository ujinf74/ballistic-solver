import argparse
import concurrent.futures
import importlib
import math
import os
import random
import statistics
import sys
import time
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
LOCAL_DEFAULT_SOLVERS = ["existing", "no_aux", "aux_fixed_point"]
ALL_SOLVERS = ["YudoTLE", "existing", "no_aux", "aux_fixed_point", "no_aux_py", "vacuum_only"]
EXPENSIVE_HIGH_ARC_SOLVERS = {"no_aux", "aux_fixed_point", "no_aux_py"}
# Benchmark policy: treat any single-case runtime above 1 second as a failure.
BENCHMARK_TIMEOUT_MS = 1000.0
PARAM_FIELDS = (
    "arcMode",
    "g",
    "wind",
    "dt",
    "tMax",
    "tolMiss",
    "beta",
    "maxIter",
    "lambdaInit",
    "lambdaMin",
    "lambdaMax",
    "lambdaUpMul",
    "lambdaDownMul",
    "lambdaTries",
    "lineSearchTries",
    "alphaMin",
    "thetaMin",
    "thetaMax",
    "broydenMinDenom",
    "missEps",
    "fdScale",
    "fdMin",
    "fdMax",
)


def bench_site_candidates():
    env_site = os.environ.get("BALLISTIC_BENCH_SITE")
    if env_site:
        yield Path(env_site)
    yield ROOT / ".bench_resume_site"
    yield ROOT / ".bench_site"


def import_bench_module(module_name, required_attrs):
    errors = []
    for site in bench_site_candidates():
        if not site.exists():
            continue
        sys.path.insert(0, str(site))
        sys.modules.pop(module_name, None)
        try:
            module = importlib.import_module(module_name)
            if all(hasattr(module, attr) for attr in required_attrs):
                return module, site
            errors.append(f"{site}: missing attrs {required_attrs}")
        except Exception as exc:
            errors.append(f"{site}: {exc}")
        finally:
            if sys.path and sys.path[0] == str(site):
                sys.path.pop(0)
            if module_name in sys.modules and not all(hasattr(sys.modules[module_name], attr) for attr in required_attrs):
                sys.modules.pop(module_name, None)
    joined = "; ".join(errors) if errors else "no candidate bench site found"
    raise RuntimeError(f"failed to import {module_name}: {joined}")


def import_workspace_ballistic_solver():
    module, site = import_bench_module("ballistic_solver", ("solve", "closest_approach", "simulate"))
    sys.path.insert(0, str(site))
    return module


def import_yudotle():
    module, site = import_bench_module("ballistic_solve", ("Environment", "Projectile", "Ballistic"))
    sys.path.insert(0, str(site))
    return module


def import_bench_variants():
    module, site = import_bench_module("bench_variants_core", ("solve_no_aux", "solve_aux_fixed_point"))
    sys.path.insert(0, str(site))
    return module


def clone_params(bs, params, *, dt=None):
    out = bs.BallisticParams()
    for field in PARAM_FIELDS:
        setattr(out, field, getattr(params, field))
    if dt is not None:
        out.dt = float(dt)
    return out


def wrap_pi(a):
    a = math.fmod(a + math.pi, 2.0 * math.pi)
    if a < 0.0:
        a += 2.0 * math.pi
    return a - math.pi


def direction(theta, phi):
    ct = math.cos(theta)
    return (ct * math.cos(phi), ct * math.sin(phi), math.sin(theta))


def norm(v):
    return math.sqrt(sum(x * x for x in v))


def vec_to_angles(r):
    n = norm(r)
    if n < 1e-12:
        return 0.0, 0.0
    dx = r[0] / n
    dy = r[1] / n
    dz = r[2] / n
    theta = math.asin(max(-1.0, min(1.0, dz)))
    phi = math.atan2(dy, dx)
    return theta, phi


def target_pos(rel_pos0, rel_vel, t):
    return tuple(rel_pos0[i] + rel_vel[i] * t for i in range(3))


def projectile_velocity_at_time(bs, case, theta, phi, t, params):
    t = max(0.0, float(t))
    full_steps = int(math.floor(t / params.dt))
    rem = t - full_steps * params.dt
    v0 = tuple(case["v0"] * x for x in direction(theta, phi))
    if full_steps > 0:
        traj = bs.simulate((0.0, 0.0, 0.0), v0, case["k_drag"], full_steps, params)
        r = tuple(float(x) for x in traj["r"][-1])
        v = tuple(float(x) for x in traj["v"][-1])
    else:
        r = (0.0, 0.0, 0.0)
        v = v0

    if rem > 1e-12:
        step_params = clone_params(bs, params, dt=rem)
        traj = bs.simulate(r, v, case["k_drag"], 1, step_params)
        v = tuple(float(x) for x in traj["v"][-1])
    return v


def initial_guess_straight_lead(case):
    t0 = norm(case["rel_pos0"]) / max(case["v0"], 1e-6)
    aim = target_pos(case["rel_pos0"], case["rel_vel"], t0)
    return vec_to_angles(aim)


def make_params(bs, arc="low", max_iter=None, lambda_tries=None, line_search_tries=None, t_max=None):
    p = bs.BallisticParams()
    p.arcMode = bs.ArcMode.High if arc == "high" else bs.ArcMode.Low
    p.g = 9.80665
    p.wind = (0.0, 0.0, 0.0)
    p.dt = 0.01
    p.tMax = 20.0 if t_max is None else float(t_max)
    p.tolMiss = 1e-2
    p.beta = 1.0
    p.maxIter = 20
    p.lambdaInit = 1e-6
    p.lambdaMin = 1e-12
    p.lambdaMax = 1e6
    p.lambdaUpMul = 10.0
    p.lambdaDownMul = 0.3
    p.lambdaTries = 6
    p.lineSearchTries = 10
    p.alphaMin = 1e-4
    p.thetaMin = -math.pi / 2.0
    p.thetaMax = math.pi / 2.0
    p.broydenMinDenom = 1e-12
    p.missEps = 1e-12
    p.fdScale = 1e-4
    p.fdMin = 1e-6
    p.fdMax = 1e-3
    if max_iter is not None:
        p.maxIter = int(max_iter)
    if lambda_tries is not None:
        p.lambdaTries = int(lambda_tries)
    if line_search_tries is not None:
        p.lineSearchTries = int(line_search_tries)
    return p


def make_cases(bs, count, params, seed, arc="low", target_mode="moving", min_hit_z=None):
    rng = random.Random(seed)
    cases = []
    attempts = 0
    max_attempts = max(count * 200, 1000)
    if arc == "high":
        theta_min = math.radians(50.0)
        theta_max = math.radians(78.0)
        t_min = 2.0
        t_max = 12.0
    else:
        theta_min = math.radians(2.0)
        theta_max = math.radians(38.0)
        t_min = 0.8
        t_max = 8.5
    while len(cases) < count:
        attempts += 1
        if attempts > max_attempts:
            raise RuntimeError(
                f"failed to generate {count} {arc}-arc cases within tMax={params.tMax:.3f}; "
                "relax angle/time ranges or increase --t-max"
            )
        speed = rng.uniform(120.0, 420.0)
        theta = rng.uniform(theta_min, theta_max)
        phi = rng.uniform(-math.pi, math.pi)
        if arc == "high":
            apex_est = speed * math.sin(theta) / max(params.g, 1e-9)
            high_t_min = max(t_min, apex_est + 0.5)
            high_t_cap = min(t_max, params.tMax)
            if high_t_min > high_t_cap:
                continue
            high_t_max = min(high_t_cap, apex_est + 8.0)
            if high_t_max < high_t_min:
                continue
            t_hit = rng.uniform(high_t_min, high_t_max)
        else:
            t_hit = rng.uniform(t_min, t_max)
        steps = max(1, round(t_hit / params.dt))
        t_hit = steps * params.dt
        rel_vel = (
            rng.uniform(-25.0, 25.0),
            rng.uniform(-25.0, 25.0),
            rng.uniform(-8.0, 8.0),
        )
        k_drag = rng.choice([0.0, 0.0005, 0.001, 0.002])
        d = direction(theta, phi)
        traj = bs.simulate((0.0, 0.0, 0.0), tuple(speed * x for x in d), k_drag, steps, params)
        hit = tuple(traj["r"][-1])
        if min_hit_z is not None and hit[2] < min_hit_z:
            continue
        if target_mode == "stationary":
            rel_vel = (0.0, 0.0, 0.0)
        rel_pos0 = tuple(hit[i] - rel_vel[i] * t_hit for i in range(3))
        cases.append(
            {
                "rel_pos0": rel_pos0,
                "rel_vel": rel_vel,
                "v0": speed,
                "k_drag": k_drag,
                "t_hit": t_hit,
                "arc": arc,
            }
        )
    return cases


def params_to_dict(params):
    return {
        "arcMode": "high" if params.arcMode == params.arcMode.High else "low",
        "g": float(params.g),
        "wind": tuple(float(x) for x in params.wind),
        "dt": float(params.dt),
        "tMax": float(params.tMax),
        "tolMiss": float(params.tolMiss),
        "beta": float(params.beta),
        "maxIter": int(params.maxIter),
        "lambdaInit": float(params.lambdaInit),
        "lambdaMin": float(params.lambdaMin),
        "lambdaMax": float(params.lambdaMax),
        "lambdaUpMul": float(params.lambdaUpMul),
        "lambdaDownMul": float(params.lambdaDownMul),
        "lambdaTries": int(params.lambdaTries),
        "lineSearchTries": int(params.lineSearchTries),
        "alphaMin": float(params.alphaMin),
        "thetaMin": float(params.thetaMin),
        "thetaMax": float(params.thetaMax),
        "broydenMinDenom": float(params.broydenMinDenom),
        "missEps": float(params.missEps),
        "fdScale": float(params.fdScale),
        "fdMin": float(params.fdMin),
        "fdMax": float(params.fdMax),
    }


def eval_candidate(bs, params, case, theta, phi):
    ca = bs.closest_approach(case["rel_pos0"], case["rel_vel"], case["v0"], case["k_drag"], theta, phi, params)
    proj_vel = projectile_velocity_at_time(bs, case, theta, phi, float(ca["tStar"]), params)
    return {
        "theta": theta,
        "phi": phi,
        "miss": float(ca["miss"]),
        "tStar": float(ca["tStar"]),
        "relMissAtStar": tuple(float(x) for x in ca["relMissAtStar"]),
        "projVzAtStar": float(proj_vel[2]),
    }


def residual_direct(bs, params, case, theta, phi):
    cand = eval_candidate(bs, params, case, theta, phi)
    aim = target_pos(case["rel_pos0"], case["rel_vel"], cand["tStar"])
    aim_corr = tuple(aim[i] - params.beta * cand["relMissAtStar"][i] for i in range(3))
    th0, ph0 = vec_to_angles(aim)
    th1, ph1 = vec_to_angles(aim_corr)
    return (th1 - th0, wrap_pi(ph1 - ph0)), cand


def jacobian_fd(bs, params, case, theta, phi, f_base):
    def pick_fd_step(x):
        h = params.fdScale * (1.0 + abs(x))
        return max(params.fdMin, min(params.fdMax, h))

    j = [[math.nan, math.nan], [math.nan, math.nan]]
    h_th = pick_fd_step(theta)
    can_minus = theta - h_th >= params.thetaMin
    can_plus = theta + h_th <= params.thetaMax

    if can_minus and can_plus:
        fp, _ = residual_direct(bs, params, case, theta + h_th, phi)
        fm, _ = residual_direct(bs, params, case, theta - h_th, phi)
        j[0][0] = (fp[0] - fm[0]) / (2.0 * h_th)
        j[1][0] = wrap_pi(fp[1] - fm[1]) / (2.0 * h_th)
    elif can_plus:
        fp, _ = residual_direct(bs, params, case, theta + h_th, phi)
        j[0][0] = (fp[0] - f_base[0]) / h_th
        j[1][0] = wrap_pi(fp[1] - f_base[1]) / h_th
    elif can_minus:
        fm, _ = residual_direct(bs, params, case, theta - h_th, phi)
        j[0][0] = (f_base[0] - fm[0]) / h_th
        j[1][0] = wrap_pi(f_base[1] - fm[1]) / h_th
    else:
        raise RuntimeError("theta finite difference step failed")

    h_ph = pick_fd_step(phi)
    fp, _ = residual_direct(bs, params, case, theta, wrap_pi(phi + h_ph))
    fm, _ = residual_direct(bs, params, case, theta, wrap_pi(phi - h_ph))
    j[0][1] = (fp[0] - fm[0]) / (2.0 * h_ph)
    j[1][1] = wrap_pi(fp[1] - fm[1]) / (2.0 * h_ph)
    return j


def solve_lm_step(j, f, lamb):
    a00 = j[0][0] * j[0][0] + j[1][0] * j[1][0] + lamb
    a01 = j[0][0] * j[0][1] + j[1][0] * j[1][1]
    a11 = j[0][1] * j[0][1] + j[1][1] * j[1][1] + lamb
    b0 = -(j[0][0] * f[0] + j[1][0] * f[1])
    b1 = -(j[0][1] * f[0] + j[1][1] * f[1])
    det = a00 * a11 - a01 * a01
    if abs(det) < 1e-18:
        return None
    dtheta = (b0 * a11 - b1 * a01) / det
    dphi = (a00 * b1 - a01 * b0) / det
    if not (math.isfinite(dtheta) and math.isfinite(dphi)):
        return None
    return dtheta, dphi


def broyden_update(j, sdx, ydf, min_denom):
    denom = sdx[0] * sdx[0] + sdx[1] * sdx[1]
    if denom <= min_denom:
        return
    jsdx0 = j[0][0] * sdx[0] + j[0][1] * sdx[1]
    jsdx1 = j[1][0] * sdx[0] + j[1][1] * sdx[1]
    u0 = (ydf[0] - jsdx0) / denom
    u1 = (ydf[1] - jsdx1) / denom
    j[0][0] += u0 * sdx[0]
    j[0][1] += u0 * sdx[1]
    j[1][0] += u1 * sdx[0]
    j[1][1] += u1 * sdx[1]


def solve_direct(bs, params, case, init_mode="vacuum"):
    if init_mode == "straight":
        theta, phi = initial_guess_straight_lead(case)
    else:
        guess = bs.vacuum_lead_initial_guess(case["rel_pos0"], case["rel_vel"], case["v0"], case["arc"], params.g)
        theta = float(guess["theta"])
        phi = float(guess["phi"])
    theta = max(params.thetaMin, min(params.thetaMax, theta))
    phi = wrap_pi(phi)
    f, cur = residual_direct(bs, params, case, theta, phi)
    best = dict(cur)
    best["iterations"] = 0
    try:
        j = jacobian_fd(bs, params, case, theta, phi, f)
    except RuntimeError:
        best["success"] = is_success(case, params, best)
        return best

    lamb = max(params.lambdaMin, min(params.lambdaMax, params.lambdaInit))
    accepted_steps = 0

    for it in range(params.maxIter):
        best["iterations"] = it + 1
        if cur["miss"] <= params.tolMiss:
            break

        accepted = False
        for _ in range(params.lambdaTries):
            step = solve_lm_step(j, f, lamb)
            if step is None:
                lamb = max(params.lambdaMin, min(params.lambdaMax, lamb * params.lambdaUpMul))
                continue

            dtheta, dphi = step
            alpha = 1.0
            chosen = None
            chosen_f = None
            miss_old = cur["miss"]

            for _ in range(params.lineSearchTries):
                th_try = max(params.thetaMin, min(params.thetaMax, theta + alpha * dtheta))
                ph_try = wrap_pi(phi + alpha * dphi)
                f_new, cand = residual_direct(bs, params, case, th_try, ph_try)
                if chosen is None or cand["miss"] < chosen["miss"]:
                    chosen = cand
                    chosen_f = f_new
                if cand["miss"] <= miss_old + params.missEps:
                    chosen = cand
                    chosen_f = f_new
                    accepted = True
                    break
                alpha *= 0.5
                if alpha < params.alphaMin:
                    break

            if not accepted and chosen is not None and chosen["miss"] < miss_old:
                accepted = True

            if not accepted or chosen is None:
                lamb = max(params.lambdaMin, min(params.lambdaMax, lamb * params.lambdaUpMul))
                continue

            accepted_steps += 1
            if chosen["miss"] < miss_old - params.missEps:
                lamb = max(params.lambdaMin, min(params.lambdaMax, lamb * params.lambdaDownMul))

            sdx = (chosen["theta"] - theta, wrap_pi(chosen["phi"] - phi))
            ydf = (chosen_f[0] - f[0], wrap_pi(chosen_f[1] - f[1]))
            broyden_update(j, sdx, ydf, params.broydenMinDenom)

            theta = chosen["theta"]
            phi = wrap_pi(chosen["phi"])
            f = chosen_f
            cur = chosen
            if cur["miss"] < best["miss"]:
                best = dict(cur)
                best["iterations"] = it + 1
            break

        if not accepted:
            break

    best["acceptedSteps"] = accepted_steps
    best["success"] = is_success(case, params, best)
    return best


def solve_vacuum_only(bs, params, case):
    guess = bs.vacuum_lead_initial_guess(case["rel_pos0"], case["rel_vel"], case["v0"], case["arc"], params.g)
    theta = float(guess["theta"])
    phi = float(guess["phi"])
    out = eval_candidate(bs, params, case, theta, phi)
    out["iterations"] = 0
    out["acceptedSteps"] = 0
    out["success"] = is_success(case, params, out)
    return out


def make_yudotle_solver(ys, k_drag, env_mode, g):
    rho = 1.225
    area = 0.0 if k_drag == 0.0 else (2.0 * k_drag / rho)
    env = ys.Environment.earth_standard()
    if env_mode == "constant":
        env = env.with_gravity(g).with_air_density(rho).with_wind_velocity(np.zeros(3))
    projectile = ys.Projectile.baseball_standard().with_mass(1.0).with_area(area).with_drag_coefficient(1.0)
    return ys.Ballistic(env, projectile)


def solve_yudotle(bs, ys, params, case, env_mode):
    solver = make_yudotle_solver(ys, case["k_drag"], env_mode, params.g)
    rel_pos0 = np.asarray(case["rel_pos0"], dtype=np.float64)
    rel_vel = np.asarray(case["rel_vel"], dtype=np.float64)

    def target_position(t):
        return rel_pos0 + rel_vel * t

    sols = []
    for method_name in ("solve_earliest", "solve_latest"):
        sol = getattr(solver, method_name)(
            target_position,
            np.zeros(3, dtype=np.float64),
            np.zeros(3, dtype=np.float64),
            float(case["v0"]),
            (0.0, params.tMax),
            0.25,
        )
        if sol is not None:
            sols.append((method_name, sol))

    if not sols:
        return {
            "theta": math.nan,
            "phi": math.nan,
            "miss": math.inf,
            "tStar": math.nan,
            "relMissAtStar": (math.nan, math.nan, math.nan),
            "projVzAtStar": math.nan,
            "iterations": math.nan,
            "acceptedSteps": math.nan,
            "success": False,
        }

    picked_out = None
    picked_meta = None
    for method_name, sol in sols:
        direction_vec = tuple(float(x) for x in sol.direction)
        theta, phi = vec_to_angles(direction_vec)
        out = eval_candidate(bs, params, case, theta, phi)
        if picked_out is None:
            picked_out = out
            picked_meta = (method_name, sol)
            continue
        if case["arc"] == "high":
            if out["theta"] > picked_out["theta"]:
                picked_out = out
                picked_meta = (method_name, sol)
        else:
            if out["theta"] < picked_out["theta"]:
                picked_out = out
                picked_meta = (method_name, sol)

    method_name, sol = picked_meta
    out = picked_out
    out["iterations"] = math.nan
    out["acceptedSteps"] = math.nan
    out["success"] = is_success(case, params, out)
    out["selectedMethod"] = method_name
    out["reportedTime"] = float(sol.time)
    out["reportedError"] = float(sol.error)
    out["reportedComputeTimeMs"] = float(sol.computation_time) * 1000.0
    return out


def solve_existing(bs, params, case):
    result = bs.solve(case["rel_pos0"], case["rel_vel"], case["v0"], case["k_drag"], params=params)
    out = eval_candidate(bs, params, case, float(result["theta"]), float(result["phi"]))
    out["iterations"] = int(result.get("iterations", 0))
    out["acceptedSteps"] = int(result.get("acceptedSteps", 0))
    out["success"] = bool(result.get("success", False)) and is_success(case, params, out)
    return out


def solve_no_aux_native(bench_variants, bs, params, case):
    result = bench_variants.solve_no_aux(
        case["rel_pos0"],
        case["rel_vel"],
        case["v0"],
        case["k_drag"],
        params_to_dict(params),
    )
    out = eval_candidate(bs, params, case, float(result["theta"]), float(result["phi"]))
    out["iterations"] = int(result.get("iterations", 0))
    out["acceptedSteps"] = int(result.get("acceptedSteps", 0))
    out["success"] = bool(result.get("success", False)) and is_success(case, params, out)
    return out


def solve_aux_fixed_point_native(bench_variants, bs, params, case):
    result = bench_variants.solve_aux_fixed_point(
        case["rel_pos0"],
        case["rel_vel"],
        case["v0"],
        case["k_drag"],
        params_to_dict(params),
    )
    out = eval_candidate(bs, params, case, float(result["theta"]), float(result["phi"]))
    out["iterations"] = int(result.get("iterations", 0))
    out["acceptedSteps"] = int(result.get("acceptedSteps", 0))
    out["success"] = bool(result.get("success", False)) and is_success(case, params, out)
    return out


def solve_case_by_name(bench_variants, bs, ys, params, case, solver_name, yudotle_env):
    if solver_name == "YudoTLE":
        return solve_yudotle(bs, ys, params, case, yudotle_env)
    if solver_name == "existing":
        return solve_existing(bs, params, case)
    if solver_name == "no_aux":
        return solve_no_aux_native(bench_variants, bs, params, case)
    if solver_name == "aux_fixed_point":
        return solve_aux_fixed_point_native(bench_variants, bs, params, case)
    if solver_name == "no_aux_py":
        return solve_direct(bs, params, case, "straight")
    if solver_name == "vacuum_only":
        return solve_vacuum_only(bs, params, case)
    raise ValueError(f"unknown solver: {solver_name}")


def is_success(case, params, out):
    if not math.isfinite(out["miss"]) or out["miss"] > params.tolMiss:
        return False
    if case["arc"] == "high":
        return math.isfinite(out["projVzAtStar"]) and out["projVzAtStar"] <= 0.0
    return True


def percentile(values, q):
    values = sorted(values)
    if not values:
        return math.nan
    return values[min(len(values) - 1, round((len(values) - 1) * q))]


def summarize(name, rows):
    times = [r["runtime_ms"] for r in rows]
    misses = [r["miss"] for r in rows]
    successes = sum(1 for r in rows if r["success"])
    finite_iters = [r["iterations"] for r in rows if isinstance(r["iterations"], (int, float)) and math.isfinite(r["iterations"])]
    finite_accepted = [r["acceptedSteps"] for r in rows if isinstance(r.get("acceptedSteps"), (int, float)) and math.isfinite(r["acceptedSteps"])]
    return {
        "name": name,
        "count": len(rows),
        "successes": successes,
        "success_rate": successes / len(rows) if rows else math.nan,
        "median_ms": statistics.median(times) if times else math.nan,
        "p95_ms": percentile(times, 0.95),
        "median_miss": statistics.median(misses) if misses else math.nan,
        "p95_miss": percentile(misses, 0.95),
        "mean_iter": statistics.mean(finite_iters) if finite_iters else math.nan,
        "median_iter": statistics.median(finite_iters) if finite_iters else math.nan,
        "mean_accepted": statistics.mean(finite_accepted) if finite_accepted else math.nan,
        "median_accepted": statistics.median(finite_accepted) if finite_accepted else math.nan,
    }


def run_case(solver_fn, case):
    t0 = time.perf_counter()
    out = solver_fn(case)
    runtime_ms = (time.perf_counter() - t0) * 1000.0
    row = dict(out)
    row["runtime_ms"] = runtime_ms
    # This is a reporting rule only; it does not stop the computation mid-flight.
    row["timeout"] = runtime_ms > BENCHMARK_TIMEOUT_MS
    if row["timeout"]:
        row["success"] = False
    return row


def run_one(name, solver_fn, cases, workers):
    t0 = time.perf_counter()
    if workers <= 1:
        rows = [run_case(solver_fn, case) for case in cases]
    else:
        with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as ex:
            rows = list(ex.map(lambda case: run_case(solver_fn, case), cases))
    wall_ms = (time.perf_counter() - t0) * 1000.0
    summary = summarize(name, rows)
    summary["wall_ms"] = wall_ms
    summary["wall_per_case_ms"] = wall_ms / len(cases) if cases else math.nan
    return summary, rows


def print_summary(summary):
    mean_iter = "nan" if not math.isfinite(summary["mean_iter"]) else f"{summary['mean_iter']:.2f}"
    median_iter = "nan" if not math.isfinite(summary["median_iter"]) else f"{summary['median_iter']:.2f}"
    mean_accepted = "nan" if not math.isfinite(summary["mean_accepted"]) else f"{summary['mean_accepted']:.2f}"
    median_accepted = "nan" if not math.isfinite(summary["median_accepted"]) else f"{summary['median_accepted']:.2f}"
    print(
        f"{summary['name']}: "
        f"success={summary['successes']}/{summary['count']} "
        f"rate={summary['success_rate'] * 100:.2f}% "
        f"median={summary['median_ms']:.3f}ms "
        f"p95={summary['p95_ms']:.3f}ms "
        f"median_miss={summary['median_miss']:.3e} "
        f"p95_miss={summary['p95_miss']:.3e} "
        f"mean_iter={mean_iter} "
        f"median_iter={median_iter} "
        f"mean_accepted={mean_accepted} "
        f"median_accepted={median_accepted} "
        f"wall={summary['wall_ms']:.3f}ms "
        f"wall_per_case={summary['wall_per_case_ms']:.3f}ms"
    )


def build_solver_table(bench_variants, bs, ys, params, yudotle_env):
    return {
        "YudoTLE": lambda case: solve_yudotle(bs, ys, params, case, yudotle_env),
        "existing": lambda case: solve_existing(bs, params, case),
        # `no_aux` is the native straight-lead direct residual solver used for fair runtime comparison.
        "no_aux": lambda case: solve_no_aux_native(bench_variants, bs, params, case),
        # Auxiliary-solution fixed-point update with alpha-only line search; no Jacobian or lambda search.
        "aux_fixed_point": lambda case: solve_aux_fixed_point_native(bench_variants, bs, params, case),
        # Keep the original Python implementation only as a benchmark/reference path.
        "no_aux_py": lambda case: solve_direct(bs, params, case, "straight"),
        "vacuum_only": lambda case: solve_vacuum_only(bs, params, case),
    }


def pick_solver_names(arc, requested):
    if requested:
        if "all" in requested:
            return ALL_SOLVERS
        return requested
    if arc == "high":
        return [name for name in LOCAL_DEFAULT_SOLVERS if name not in EXPENSIVE_HIGH_ARC_SOLVERS]
    return LOCAL_DEFAULT_SOLVERS


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--count", type=int, default=2000)
    parser.add_argument("--seed", type=int, default=20260503)
    parser.add_argument("--arc", choices=["low", "high"], default="low")
    parser.add_argument("--yudotle-env", choices=["earth", "constant"], default="constant")
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--t-max", type=float, default=None)
    parser.add_argument("--target-mode", choices=["moving", "stationary"], default="moving")
    parser.add_argument("--min-hit-z", type=float, default=None)
    parser.add_argument(
        "--solvers",
        nargs="+",
        choices=["all", "YudoTLE", "existing", "no_aux", "aux_fixed_point", "no_aux_py", "vacuum_only"],
        default=None,
    )
    parser.add_argument("--max-iter", type=int, default=None)
    parser.add_argument("--lambda-tries", type=int, default=None)
    parser.add_argument("--line-search-tries", type=int, default=None)
    args = parser.parse_args()

    bs = import_workspace_ballistic_solver()
    bench_variants = import_bench_variants()
    params = make_params(bs, args.arc, args.max_iter, args.lambda_tries, args.line_search_tries, args.t_max)
    cases = make_cases(
        bs,
        args.count,
        params,
        args.seed,
        arc=args.arc,
        target_mode=args.target_mode,
        min_hit_z=args.min_hit_z,
    )
    solver_names = pick_solver_names(args.arc, args.solvers)
    ys = import_yudotle() if "YudoTLE" in solver_names else None
    solver_table = build_solver_table(bench_variants, bs, ys, params, args.yudotle_env)

    workers = int(args.workers)
    if workers == 0:
        workers = max(1, os.cpu_count() or 1)
    workers = max(1, workers)

    if args.arc == "high" and args.solvers is None:
        skipped = ", ".join(sorted(EXPENSIVE_HIGH_ARC_SOLVERS))
        print(f"info: high-arc default skips expensive solvers: {skipped}. Use --solvers all or explicit names to include them.")
    elif args.arc == "high" and any(name in EXPENSIVE_HIGH_ARC_SOLVERS for name in solver_names):
        chosen = ", ".join(name for name in solver_names if name in EXPENSIVE_HIGH_ARC_SOLVERS)
        print(f"info: running expensive high-arc solvers: {chosen}")

    for name in solver_names:
        fn = solver_table[name]
        summary, _ = run_one(name, fn, cases, workers)
        print_summary(summary)


if __name__ == "__main__":
    main()
