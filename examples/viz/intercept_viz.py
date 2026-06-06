"""Visualizations for the moving-target intercept evaluation.

Produces (into examples/viz/renders/):
  * intercept_miss_vs_range.png -- real miss distance vs closing range for the
    tracker, the ideal (true-state) solution, and a naive CV baseline, across the
    robotics / security / agriculture / defense scenarios.
  * intercept_engagement.gif    -- top-down animation of a clean straight-inbound
    intercept (true path, noisy track, predicted intercept, projectile tracers,
    hits). The hard weaving (defense) case is the limit shown in the curves.

    python examples/viz/intercept_viz.py

Reuses the scenarios and real-miss simulation from benchmarks/intercept_eval.
matplotlib is an optional demo-only dependency.
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
            rows.append((
                R,
                eng.real_miss(v0, k_drag, trk["theta"], trk["phi"], true_at),
                eng.real_miss(v0, k_drag, ideal["theta"], ideal["phi"], true_at),
                cv_miss,
            ))
        prev, prev_t = meas, t
    return rows


def plot_miss_vs_range():
    fig, axes = plt.subplots(2, 2, figsize=(13, 9), sharey=True)
    for ax, scn_fn in zip(axes.flat, eng.SCENARIOS):
        scn = scn_fn()
        rows = burst_curve(scn)
        if not rows:
            ax.set_title(f"{scn['name']} ({scn['domain']}) - no data")
            continue
        rng = [r[0] for r in rows]
        ax.semilogy(rng, [r[2] for r in rows], label="ideal (true state)", color="tab:green", lw=1.6)
        ax.semilogy(rng, [r[1] for r in rows], label="tracker (noisy track)", color="tab:blue", lw=1.6)
        ax.semilogy(rng, [r[3] for r in rows], label="naive CV", color="tab:red", lw=1.2, alpha=0.8)
        ax.axhspan(0.0, 1.0, color="gray", alpha=0.18)
        ax.set_title(f"{scn['name']}  ({scn['domain']}, v0={scn['v0']:.0f} m/s)")
        ax.set_xlabel("closing range (m)")
        ax.invert_xaxis()
        ax.grid(True, which="both", alpha=0.3)
    axes[0, 0].set_ylabel("real miss (m)")
    axes[1, 0].set_ylabel("real miss (m)")
    axes[0, 0].legend(loc="upper right", fontsize=8)
    fig.suptitle("Moving-target intercept: real miss vs closing range (shaded = within 1 m)")
    fig.tight_layout()
    path = OUT / "intercept_miss_vs_range.png"
    fig.savefig(path, dpi=120)
    plt.close(fig)
    return path


# --- clean headline engagement: a straight inbound target (no maneuver) ---
def demo_scenario():
    s, v = (300.0, 40.0, 30.0), (-45.0, -6.0, 0.0)
    return {
        "name": "straight inbound", "v0": 160.0, "k_drag": 0.006,
        "pos": lambda t: (s[0] + v[0] * t, s[1] + v[1] * t, s[2]),
        "fire_range": 240.0, "min_range": 35.0, "tmax": 5.6,
    }


def projectile_xy(v0, k_drag, theta, phi, dt=0.01, tmax=4.0):
    d = eng.direction(theta, phi)
    vv = (v0 * d[0], v0 * d[1], v0 * d[2])
    P = bs.params_preset("balanced")
    P.dt = dt
    traj = bs.simulate((0.0, 0.0, 0.0), vv, k_drag, int(tmax / dt), P)
    return traj["r"]


def animate_engagement(noise=1.0, q=3.0, seed=11, hit_thresh=2.0):
    scn = demo_scenario()
    v0, k_drag = scn["v0"], scn["k_drag"]
    rng = random.Random(seed)
    tracker = bs.TargetTracker(processNoise=q, measNoise=max(noise * noise, 1e-6))
    P = eng.solve_params()

    frame_dt = 0.04
    fire_period = 0.16
    times = [k * frame_dt for k in range(int(scn["tmax"] / frame_dt))]
    true_pts = [scn["pos"](t) for t in times]
    meas_pts = [tuple(p[i] + rng.gauss(0.0, noise) for i in range(3)) for p in true_pts]

    shots = []  # (launch_t, path_r, dt, flight, hit, pred_xy)
    last_fire = -1e9
    for k, t in enumerate(times):
        tracker.update(t, meas_pts[k])
        R = math.dist(true_pts[k], (0.0, 0.0, 0.0))
        if R < scn["fire_range"] and R > scn["min_range"] and (t - last_fire) >= fire_period:
            last_fire = t
            sol = tracker.solve(v0, k_drag, params=P)
            tof = sol["tStar"]
            pr = projectile_xy(v0, k_drag, sol["theta"], sol["phi"], dt=0.01, tmax=tof)
            true_at = lambda s, t=t: scn["pos"](t + s)
            miss = eng.real_miss(v0, k_drag, sol["theta"], sol["phi"], true_at)
            pi = tracker.predict(tof)
            shots.append((t, pr, 0.01, tof, miss <= hit_thresh, (pi[0], pi[1])))

    fig, ax = plt.subplots(figsize=(8, 6))
    ax.set_aspect("equal")
    ax.set_xlim(-40, 340)
    ax.set_ylim(-120, 120)
    ax.set_xlabel("down-range x (m)")
    ax.set_ylabel("cross-range y (m)")
    ax.grid(True, alpha=0.3)
    ax.plot([p[0] for p in true_pts], [p[1] for p in true_pts], color="tab:green", lw=1.0, alpha=0.4, label="true path")
    ax.plot(0, 0, marker="^", color="black", ms=12, label="launcher")
    tgt_dot, = ax.plot([], [], "o", color="tab:green", ms=10, label="target")
    meas_sc = ax.scatter([], [], s=12, color="tab:red", alpha=0.5, label="noisy track")
    pred_dot, = ax.plot([], [], "x", color="tab:orange", ms=11, mew=2, label="predicted intercept")
    tracer_lines = [ax.plot([], [], color="dimgray", lw=0.8, alpha=0.6)[0] for _ in range(8)]
    proj_dot, = ax.plot([], [], "o", color="tab:blue", ms=4, label="rounds in flight")
    hits_x, hits_y = [], []
    hits_sc = ax.scatter([], [], marker="*", s=160, color="lime", edgecolor="darkgreen", label="hit (<%.0f m)" % hit_thresh, zorder=5)
    txt = ax.text(0.02, 0.97, "", transform=ax.transAxes, va="top", fontsize=9, family="monospace")
    ax.legend(loc="lower right", fontsize=8)

    n_hits = [0]
    scored = set()

    def update(frame):
        t = times[frame]
        ax.set_title("straight inbound - clean lead intercept")
        tgt_dot.set_data([true_pts[frame][0]], [true_pts[frame][1]])
        lo = max(0, frame - 18)
        meas_sc.set_offsets([[meas_pts[i][0], meas_pts[i][1]] for i in range(lo, frame + 1)])

        active = [(tf, pr, d, fl, hit, pxy) for (tf, pr, d, fl, hit, pxy) in shots if 0.0 <= t - tf <= fl]
        if active:
            pred_dot.set_data([active[-1][5][0]], [active[-1][5][1]])

        px, py = [], []
        for li in tracer_lines:
            li.set_data([], [])
        for j, (tf, pr, d, fl, hit, pxy) in enumerate(active[-len(tracer_lines):]):
            s = t - tf
            idx = min(len(pr) - 1, int(s / d))
            px.append(pr[idx][0])
            py.append(pr[idx][1])
            tracer_lines[j].set_data([pr[i][0] for i in range(idx + 1)], [pr[i][1] for i in range(idx + 1)])
        proj_dot.set_data(px, py)

        for si, (tf, pr, d, fl, hit, pxy) in enumerate(shots):
            if hit and si not in scored and t - tf >= fl:
                scored.add(si)
                n_hits[0] += 1
                hits_x.append(true_pts[min(frame, len(true_pts) - 1)][0])
                hits_y.append(true_pts[min(frame, len(true_pts) - 1)][1])
        if hits_x:
            hits_sc.set_offsets(list(zip(hits_x, hits_y)))

        R = math.dist(true_pts[frame], (0.0, 0.0, 0.0))
        fired = sum(1 for sh in shots if sh[0] <= t)
        pk = (100.0 * n_hits[0] / fired) if fired else 0.0
        txt.set_text(f"t={t:4.1f}s  range={R:5.0f}m  fired={fired}  hits<{hit_thresh:.0f}m={n_hits[0]} ({pk:.0f}%)")
        return tgt_dot, meas_sc, pred_dot, proj_dot, hits_sc, txt

    anim = FuncAnimation(fig, update, frames=len(times), blit=False, interval=40)
    path = OUT / "intercept_engagement.gif"
    anim.save(path, writer=PillowWriter(fps=25))
    plt.close(fig)
    return path, len(shots), n_hits[0]


def main():
    p1 = plot_miss_vs_range()
    print(f"wrote {p1}")
    p2, fired, hits = animate_engagement()
    print(f"wrote {p2}  (fired={fired}, hits={hits})")


if __name__ == "__main__":
    main()
