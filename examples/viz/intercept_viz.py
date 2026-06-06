"""Visualizations for the moving-target intercept evaluation.

Produces (into examples/viz/renders/):
  * intercept_miss_vs_range.png -- real miss vs closing range (binned median +
    inter-quartile band) for the tracker, ideal (true-state), and naive CV
    baseline, across the robotics / security / agriculture / defense scenarios.
  * intercept_engagement.gif    -- a tracking-display animation of a clean
    straight-inbound intercept: noisy measurements, the tracker's denoised
    estimate, the predicted-intercept lead line, projectile tracers, and hits.

    python examples/viz/intercept_viz.py

Reuses scenarios and the real-miss simulation from benchmarks/intercept_eval.
matplotlib is an optional demo-only dependency.
"""

import math
import random
import sys
from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

import ballistic_solver as bs
import benchmarks.intercept_eval as eng

OUT = Path(__file__).resolve().parent / "renders"
OUT.mkdir(exist_ok=True)


def burst_curve(scn, *, noise=1.0, q=5.0, seed=11):
    rng = random.Random(seed)
    tracker = bs.TargetTracker(processNoise=q, measNoise=max(noise * noise, 1e-6))
    P = eng.solve_params()
    v0, k_drag = scn["v0"], scn["k_drag"]
    min_range = max(20.0, 0.05 * scn["fire_range"])
    dt = 0.01
    rows = []
    prev = prev_t = None
    for k in range(int(scn["tmax"] / dt)):
        t = k * dt
        true = scn["pos"](t)
        meas = tuple(true[i] + rng.gauss(0.0, noise) for i in range(3))
        tracker.update(t, meas)
        R = math.dist(true, (0.0, 0.0, 0.0))
        if R < scn["fire_range"] and R > min_range:
            true_at = lambda s, t=t: scn["pos"](t + s)
            trk = tracker.solve(v0, k_drag, params=P)
            ideal = bs.solve_accel(scn["pos"](t), scn["vel"](t), scn["acc"](t), v0, k_drag, params=P)
            cv_miss = math.nan
            if prev is not None:
                dtm = t - prev_t
                vcv = tuple((meas[i] - prev[i]) / dtm for i in range(3))
                cv = bs.solve(meas, vcv, v0, k_drag, params=P)
                cv_miss = eng.real_miss(v0, k_drag, cv["theta"], cv["phi"], true_at)
            rows.append((R, eng.real_miss(v0, k_drag, trk["theta"], trk["phi"], true_at),
                         eng.real_miss(v0, k_drag, ideal["theta"], ideal["phi"], true_at), cv_miss))
        prev, prev_t = meas, t
    return rows


def _binned(xs, ys, nb=16):
    xs, ys = np.asarray(xs, float), np.asarray(ys, float)
    good = np.isfinite(xs) & np.isfinite(ys)
    xs, ys = xs[good], ys[good]
    if xs.size == 0:
        return np.array([]), np.array([]), np.array([]), np.array([])
    edges = np.linspace(xs.min(), xs.max(), nb + 1)
    cx, med, lo, hi = [], [], [], []
    for i in range(nb):
        m = (xs >= edges[i]) & (xs <= edges[i + 1] if i == nb - 1 else xs < edges[i + 1])
        if not m.any():
            continue
        v = ys[m]
        cx.append(0.5 * (edges[i] + edges[i + 1]))
        med.append(np.median(v)); lo.append(np.percentile(v, 25)); hi.append(np.percentile(v, 75))
    return map(np.asarray, (cx, med, lo, hi))


def plot_miss_vs_range():
    fig, axes = plt.subplots(1, 4, figsize=(18, 4.6), sharey=True)
    for ax, scn_fn in zip(axes, eng.SCENARIOS):
        scn = scn_fn()
        rows = burst_curve(scn)
        if not rows:
            ax.set_title(f"{scn['name']} ({scn['domain']})")
            continue
        rng = [r[0] for r in rows]
        cx, med_i, _, _ = _binned(rng, [r[2] for r in rows])
        ax.plot(cx, med_i, color="tab:green", lw=2.0, label="ideal (true state)")
        cx, med_t, lo_t, hi_t = _binned(rng, [r[1] for r in rows])
        ax.fill_between(cx, lo_t, hi_t, color="tab:blue", alpha=0.20)
        ax.plot(cx, med_t, color="tab:blue", lw=2.2, label="tracker (noisy track)")
        cx, med_c, _, _ = _binned(rng, [r[3] for r in rows])
        ax.plot(cx, med_c, color="tab:red", lw=1.8, ls="--", label="naive CV")
        ax.axhspan(0.0, 1.0, color="0.6", alpha=0.25)
        ax.set_yscale("log")
        ax.set_title(f"{scn['name']}\n({scn['domain']}, v0={scn['v0']:.0f} m/s)", fontsize=10)
        ax.set_xlabel("closing range (m)")
        ax.invert_xaxis()
        ax.grid(True, which="both", alpha=0.25)
    axes[0].set_ylabel("real miss (m)  [log]")
    axes[0].legend(loc="upper right", fontsize=8, framealpha=0.9)
    fig.suptitle("Moving-target intercept: real miss vs closing range  (binned median + IQR; shaded floor = within 1 m)", fontsize=12)
    fig.tight_layout()
    path = OUT / "intercept_miss_vs_range.png"
    fig.savefig(path, dpi=120)
    plt.close(fig)
    return path


def demo_scenario():
    s, v = (300.0, 40.0, 30.0), (-45.0, -6.0, 0.0)
    return {
        "name": "straight inbound", "v0": 160.0, "k_drag": 0.006,
        "pos": lambda t: (s[0] + v[0] * t, s[1] + v[1] * t, s[2]),
        "fire_range": 240.0, "min_range": 30.0, "tmax": 5.8,
    }


def projectile_xy(v0, k_drag, theta, phi, dt=0.01, tmax=4.0):
    d = eng.direction(theta, phi)
    P = bs.params_preset("balanced")
    P.dt = dt
    traj = bs.simulate((0.0, 0.0, 0.0), (v0 * d[0], v0 * d[1], v0 * d[2]), k_drag, int(tmax / dt), P)
    return traj["r"]


def animate_engagement(noise=1.0, q=3.0, seed=11, hit_thresh=2.0):
    scn = demo_scenario()
    v0, k_drag = scn["v0"], scn["k_drag"]
    rng = random.Random(seed)
    tracker = bs.TargetTracker(processNoise=q, measNoise=max(noise * noise, 1e-6))
    P = eng.solve_params()

    frame_dt, fire_period = 0.045, 0.16
    times = [k * frame_dt for k in range(int(scn["tmax"] / frame_dt))]
    true_pts = [scn["pos"](t) for t in times]
    meas_pts, est_pts = [], []
    shots, pulses = [], []  # shot: (tf, path, dt, flight, hit, pred_xy); pulse: (frame, x, y)
    last_fire = -1e9
    for k, t in enumerate(times):
        m = tuple(true_pts[k][i] + rng.gauss(0.0, noise) for i in range(3))
        tracker.update(t, m)
        meas_pts.append(m)
        est_pts.append(tuple(tracker.position))
        R = math.dist(true_pts[k], (0.0, 0.0, 0.0))
        if R < scn["fire_range"] and R > scn["min_range"] and (t - last_fire) >= fire_period:
            last_fire = t
            sol = tracker.solve(v0, k_drag, params=P)
            tof = sol["tStar"]
            pr = projectile_xy(v0, k_drag, sol["theta"], sol["phi"], dt=0.01, tmax=tof)
            true_at = lambda s, t=t: scn["pos"](t + s)
            hit = eng.real_miss(v0, k_drag, sol["theta"], sol["phi"], true_at) <= hit_thresh
            pi = tracker.predict(tof)
            shots.append((t, pr, 0.01, tof, hit, (pi[0], pi[1]), k))

    plt.style.use("dark_background")
    fig, ax = plt.subplots(figsize=(9, 6.2))
    ax.set_aspect("equal")
    ax.set_xlim(-30, 330); ax.set_ylim(-90, 90)
    ax.set_xlabel("down-range x (m)"); ax.set_ylabel("cross-range y (m)")
    ax.grid(True, alpha=0.15)
    ax.plot([p[0] for p in true_pts], [p[1] for p in true_pts], color="0.4", lw=1.0, alpha=0.6)
    ax.plot(0, 0, marker="^", color="white", ms=13, label="launcher")
    trail, = ax.plot([], [], color="limegreen", lw=2.0, alpha=0.5)
    tgt_dot, = ax.plot([], [], "o", color="limegreen", ms=11, label="target (true)")
    meas_sc = ax.scatter([], [], s=14, color="tomato", alpha=0.45, label="noisy measurements")
    est_dot, = ax.plot([], [], "o", color="deepskyblue", ms=7, label="tracker estimate")
    lead_line, = ax.plot([], [], color="gold", lw=1.3, ls="--", alpha=0.9, label="lead to predicted intercept")
    pred_dot, = ax.plot([], [], "x", color="gold", ms=11, mew=2)
    tracers = [ax.plot([], [], color="0.7", lw=0.7, alpha=0.55)[0] for _ in range(7)]
    proj_dot, = ax.plot([], [], "o", color="white", ms=4, label="rounds in flight")
    pulse_sc = ax.scatter([], [], s=[], facecolors="none", edgecolors="lime", linewidths=2.0)
    hits_sc = ax.scatter([], [], marker="*", s=150, color="lime", edgecolor="white", zorder=6, label="hit (<%.0f m)" % hit_thresh)
    txt = ax.text(0.015, 0.97, "", transform=ax.transAxes, va="top", fontsize=10, family="monospace", color="white")
    ax.legend(loc="lower right", fontsize=8, framealpha=0.3)

    hits_xy, scored, n_hits = [], set(), [0]

    def update(frame):
        t = times[frame]
        ax.set_title("position-only tracking + lead intercept", color="white")
        tgt_dot.set_data([true_pts[frame][0]], [true_pts[frame][1]])
        tl = range(max(0, frame - 14), frame + 1)
        trail.set_data([true_pts[i][0] for i in tl], [true_pts[i][1] for i in tl])
        lo = max(0, frame - 14)
        meas_sc.set_offsets([[meas_pts[i][0], meas_pts[i][1]] for i in range(lo, frame + 1)])
        est_dot.set_data([est_pts[frame][0]], [est_pts[frame][1]])

        active = [s for s in shots if 0.0 <= t - s[0] <= s[3]]
        if active:
            pxy = active[-1][5]
            pred_dot.set_data([pxy[0]], [pxy[1]])
            lead_line.set_data([0.0, pxy[0]], [0.0, pxy[1]])
        px, py = [], []
        for li in tracers:
            li.set_data([], [])
        for j, s in enumerate(active[-len(tracers):]):
            idx = min(len(s[1]) - 1, int((t - s[0]) / s[2]))
            px.append(s[1][idx][0]); py.append(s[1][idx][1])
            tracers[j].set_data([s[1][i][0] for i in range(idx + 1)], [s[1][i][1] for i in range(idx + 1)])
        proj_dot.set_data(px, py)

        for si, s in enumerate(shots):
            if s[4] and si not in scored and t - s[0] >= s[3]:
                scored.add(si); n_hits[0] += 1
                hx, hy = true_pts[min(frame, len(true_pts) - 1)][0], true_pts[min(frame, len(true_pts) - 1)][1]
                hits_xy.append((hx, hy)); pulses.append((frame, hx, hy))
        if hits_xy:
            hits_sc.set_offsets(hits_xy)
        live = [(x, y, frame - f0) for (f0, x, y) in pulses if 0 <= frame - f0 <= 6]
        if live:
            pulse_sc.set_offsets([[x, y] for (x, y, a) in live])
            pulse_sc.set_sizes([120 + 260 * a for (x, y, a) in live])
        else:
            pulse_sc.set_offsets(np.empty((0, 2)))
            pulse_sc.set_sizes([])

        R = math.dist(true_pts[frame], (0.0, 0.0, 0.0))
        fired = sum(1 for s in shots if s[0] <= t)
        pk = (100.0 * n_hits[0] / fired) if fired else 0.0
        txt.set_text(f"t={t:4.1f}s   range={R:5.0f} m   fired={fired:2d}   hits<{hit_thresh:.0f}m={n_hits[0]:2d} ({pk:3.0f}%)")
        return ()

    anim = FuncAnimation(fig, update, frames=len(times), blit=False, interval=45)
    path = OUT / "intercept_engagement.gif"
    anim.save(path, writer=PillowWriter(fps=22))
    plt.close(fig)
    plt.style.use("default")
    return path, len(shots), n_hits[0]


def main():
    p1 = plot_miss_vs_range()
    print(f"wrote {p1}")
    p2, fired, hits = animate_engagement()
    print(f"wrote {p2}  (fired={fired}, hits={hits})")


if __name__ == "__main__":
    main()
