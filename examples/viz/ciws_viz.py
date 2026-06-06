"""Visualizations for the CIWS / gun-AA engagement evaluation.

Produces (into examples/viz/renders/):
  * ciws_miss_vs_range.png  -- real miss distance vs closing range for the
    tracker, the ideal (true-state) solution, and a naive CV baseline.
  * ciws_engagement.gif     -- top-down animation of a weaving-target engagement
    (true path, noisy radar, predicted intercept, projectile tracers, hits).

    python examples/viz/ciws_viz.py

Reuses the scenarios and the real-miss simulation from benchmarks/ciws_eval.
matplotlib is required only for these demos, not for the library itself.
"""

import math
import random
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

import ballistic_solver as bs
import benchmarks.ciws_eval as eng

OUT = Path(__file__).resolve().parent / "renders"
OUT.mkdir(exist_ok=True)


def burst_curve(scn, *, tmax, noise=1.0, q=5.0, seed=11):
    rng = random.Random(seed)
    tracker = bs.TargetTracker(processNoise=q, measNoise=max(noise * noise, 1e-6))
    P = eng.solve_params()
    dt = 0.01
    rows = []
    prev = prev_t = None
    for k in range(int(tmax / dt)):
        t = k * dt
        true = scn["pos"](t)
        meas = tuple(true[i] + rng.gauss(0.0, noise) for i in range(3))
        tracker.update(t, meas)
        R = math.dist(true, (0.0, 0.0, 0.0))
        if R < scn["fire_range"] and R > 120.0:
            true_at = lambda s, t=t: scn["pos"](t + s)
            trk = tracker.solve(eng.V0, eng.K_DRAG, params=P)
            ideal = bs.solve_accel(scn["pos"](t), scn["vel"](t), scn["acc"](t), eng.V0, eng.K_DRAG, params=P)
            cv_miss = math.nan
            if prev is not None:
                dtm = t - prev_t
                vcv = tuple((meas[i] - prev[i]) / dtm for i in range(3))
                cv = bs.solve(meas, vcv, eng.V0, eng.K_DRAG, params=P)
                cv_miss = eng.real_miss(cv["theta"], cv["phi"], true_at)
            rows.append((
                R,
                eng.real_miss(trk["theta"], trk["phi"], true_at),
                eng.real_miss(ideal["theta"], ideal["phi"], true_at),
                cv_miss,
            ))
        prev, prev_t = meas, t
    return rows


def plot_miss_vs_range():
    scenarios = [
        ("sea-skimmer, straight", eng.straight_inbound(), 9.0),
        ("sea-skimmer, weaving ~16g", eng.weaving_inbound(), 9.0),
        ("aircraft, crossing", eng.crossing(), 8.0),
    ]
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5), sharey=True)
    for ax, (name, scn, tmax) in zip(axes, scenarios):
        rows = burst_curve(scn, tmax=tmax)
        rng = [r[0] for r in rows]
        ax.semilogy(rng, [r[2] for r in rows], label="ideal (true state)", color="tab:green", lw=1.6)
        ax.semilogy(rng, [r[1] for r in rows], label="tracker (radar)", color="tab:blue", lw=1.6)
        ax.semilogy(rng, [r[3] for r in rows], label="naive CV", color="tab:red", lw=1.2, alpha=0.8)
        ax.axhspan(0.0, 5.0, color="gray", alpha=0.15)
        ax.set_title(name)
        ax.set_xlabel("closing range (m)")
        ax.invert_xaxis()
        ax.grid(True, which="both", alpha=0.3)
    axes[0].set_ylabel("real miss distance (m)")
    axes[0].legend(loc="upper right", fontsize=8)
    fig.suptitle("CIWS engagement: real miss vs closing range (shaded = within 5 m)")
    fig.tight_layout()
    path = OUT / "ciws_miss_vs_range.png"
    fig.savefig(path, dpi=120)
    plt.close(fig)
    return path


def projectile_xy(theta, phi, dt=0.01, tmax=4.0):
    d = eng.direction(theta, phi)
    vv = (eng.V0 * d[0], eng.V0 * d[1], eng.V0 * d[2])
    P = bs.params_preset("balanced")
    P.dt = dt
    traj = bs.simulate((0.0, 0.0, 0.0), vv, eng.K_DRAG, int(tmax / dt), P)
    return traj["r"], traj["t"]


def animate_engagement(noise=1.0, q=5.0, seed=11, hit_thresh=5.0):
    scn = eng.weaving_inbound()
    rng = random.Random(seed)
    tracker = bs.TargetTracker(processNoise=q, measNoise=max(noise * noise, 1e-6))
    P = eng.solve_params()

    frame_dt = 0.04
    tmax = 8.5
    fire_period = 0.18

    times = [k * frame_dt for k in range(int(tmax / frame_dt))]
    true_pts = [scn["pos"](t) for t in times]
    meas_pts = [tuple(p[i] + rng.gauss(0.0, noise) for i in range(3)) for p in true_pts]

    shots = []  # (launch_t, path_r, path_t, hit)
    last_fire = -1e9
    for k, t in enumerate(times):
        tracker.update(t, meas_pts[k])
        R = math.dist(true_pts[k], (0.0, 0.0, 0.0))
        if R < scn["fire_range"] and R > 120.0 and (t - last_fire) >= fire_period:
            last_fire = t
            sol = tracker.solve(eng.V0, eng.K_DRAG, params=P)
            pr, pt = projectile_xy(sol["theta"], sol["phi"])
            true_at = lambda s, t=t: scn["pos"](t + s)
            miss = eng.real_miss(sol["theta"], sol["phi"], true_at)
            pi = tracker.predict(sol["tStar"])
            shots.append((t, pr, pt, miss <= hit_thresh, (pi[0], pi[1])))

    fig, ax = plt.subplots(figsize=(8, 6))
    ax.set_aspect("equal")
    ax.set_xlim(-200, 2800)
    ax.set_ylim(-400, 600)
    ax.set_xlabel("down-range x (m)")
    ax.set_ylabel("cross-range y (m)")
    ax.grid(True, alpha=0.3)
    ax.plot([p[0] for p in true_pts], [p[1] for p in true_pts], color="tab:green", lw=1.0, alpha=0.5, label="true path")
    ax.plot(0, 0, marker="^", color="black", ms=12, label="gun")
    tgt_dot, = ax.plot([], [], "o", color="tab:green", ms=9, label="target")
    meas_sc = ax.scatter([], [], s=10, color="tab:red", alpha=0.5, label="radar meas")
    pred_dot, = ax.plot([], [], "x", color="tab:blue", ms=10, mew=2, label="predicted intercept")
    proj_dot, = ax.plot([], [], ".", color="dimgray", ms=5, label="rounds in flight")
    flash_dot, = ax.plot([], [], "o", color="lime", ms=14, alpha=0.0)
    hit_marks_x, hit_marks_y = [], []
    hits_sc = ax.scatter([], [], marker="*", s=120, color="lime", label="hit (<5 m)")
    txt = ax.text(0.02, 0.97, "", transform=ax.transAxes, va="top", fontsize=9, family="monospace")
    ax.legend(loc="lower left", fontsize=8)

    n_hits = [0]

    def update(frame):
        t = times[frame]
        ax.set_title("weaving target (~16g) - continuous fire")
        tgt_dot.set_data([true_pts[frame][0]], [true_pts[frame][1]])
        lo = max(0, frame - 25)
        meas_sc.set_offsets([[meas_pts[i][0], meas_pts[i][1]] for i in range(lo, frame + 1)])

        latest = [sh for sh in shots if sh[0] <= t]
        if latest:
            pxy = latest[-1][4]
            pred_dot.set_data([pxy[0]], [pxy[1]])

        # predicted intercept from current tracker (recompute lightly)
        # use the most recent shot's solution proxy: predict at its tof
        proj_x, proj_y = [], []
        for (tf, pr, pt, hit, pxy) in shots:
            s = t - tf
            if 0.0 <= s <= pt[-1]:
                idx = min(len(pt) - 1, int(s / (pt[1] - pt[0])))
                proj_x.append(pr[idx][0])
                proj_y.append(pr[idx][1])
                if hit and abs(s - pt[idx]) < frame_dt and (pr[idx][0], pr[idx][1]) not in zip(hit_marks_x, hit_marks_y):
                    pass
        proj_dot.set_data(proj_x, proj_y)

        # mark hits at the frame a hitting round reaches its target time (approx at end of flight near target)
        for (tf, pr, pt, hit, pxy) in shots:
            s = t - tf
            if hit and 0.0 <= s <= pt[-1]:
                # closest approach roughly where projectile nearest target now
                tp = scn["pos"](t)
                idx = min(len(pt) - 1, int(s / (pt[1] - pt[0])))
                if math.dist((pr[idx][0], pr[idx][1], pr[idx][2]), tp) <= 8.0:
                    if (round(tp[0]), round(tp[1])) not in set(zip([round(x) for x in hit_marks_x], [round(y) for y in hit_marks_y])):
                        hit_marks_x.append(tp[0])
                        hit_marks_y.append(tp[1])
                        n_hits[0] += 1
        if hit_marks_x:
            hits_sc.set_offsets(list(zip(hit_marks_x, hit_marks_y)))

        R = math.dist(true_pts[frame], (0.0, 0.0, 0.0))
        txt.set_text(f"t={t:4.1f}s  range={R:6.0f}m  fired={sum(1 for sh in shots if sh[0] <= t)}  hits<5m={n_hits[0]}")
        return tgt_dot, meas_sc, pred_dot, proj_dot, hits_sc, txt

    anim = FuncAnimation(fig, update, frames=len(times), blit=False, interval=40)
    path = OUT / "ciws_engagement.gif"
    anim.save(path, writer=PillowWriter(fps=25))
    plt.close(fig)
    return path, len(shots), n_hits[0]


def main():
    p1 = plot_miss_vs_range()
    print(f"wrote {p1}")
    p2, fired, hits = animate_engagement()
    print(f"wrote {p2}  (fired={fired}, hits<5m={hits})")


if __name__ == "__main__":
    main()
