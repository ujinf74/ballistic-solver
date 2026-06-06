"""Lead-prediction accuracy: CA-Kalman tracker vs. a linear (CV) baseline.

For each synthetic track we feed noisy position-only measurements to
`bs.TargetTracker` and, at every step, predict the target position `tau` seconds
ahead. The error against ground truth is compared with a constant-velocity
baseline that estimates velocity from the last two noisy measurements -- the
"velocity linear interpolation" the tracker is meant to beat.

    python benchmarks/predictor_eval.py

Lower RMS is better. The tracker should match the CV baseline on a
constant-velocity track and beat it on accelerating / turning (maneuvering)
tracks, where the second-order state estimate captures curvature the CV
baseline cannot.
"""

import math
import random


import ballistic_solver as bs


def rms(errs):
    return math.sqrt(sum(e * e for e in errs) / len(errs)) if errs else float("nan")


def evaluate(track, *, dt=0.05, duration=8.0, noise=0.5, tau=1.0, q=20.0, seed=7):
    rng = random.Random(seed)
    tracker = bs.TargetTracker(processNoise=q, measNoise=noise * noise)

    kf_err = []
    cv_err = []
    prev_meas = None
    prev_t = None

    steps = int(duration / dt) + 1
    for k in range(steps):
        t = k * dt
        true = track(t)
        meas = tuple(true[i] + rng.gauss(0.0, noise) for i in range(3))
        tracker.update(t, meas)

        future = track(t + tau)
        kf_err.append(math.dist(tracker.predict(tau), future))

        if prev_meas is not None:
            dtm = t - prev_t
            v_cv = tuple((meas[i] - prev_meas[i]) / dtm for i in range(3))
            p_cv = tuple(meas[i] + v_cv[i] * tau for i in range(3))
            cv_err.append(math.dist(p_cv, future))

        prev_meas = meas
        prev_t = t

    warm = int(1.0 / dt)  # let the filter converge before scoring
    return rms(kf_err[warm:]), rms(cv_err[warm:])


def make_cv():
    p0, v = (140.0, 40.0, 8.0), (6.0, -4.0, 1.0)
    return lambda t: (p0[0] + v[0] * t, p0[1] + v[1] * t, p0[2] + v[2] * t)


def make_ca():
    p0, v, a = (140.0, 40.0, 8.0), (6.0, -4.0, 1.0), (0.8, 0.6, -0.2)
    return lambda t: (
        p0[0] + v[0] * t + 0.5 * a[0] * t * t,
        p0[1] + v[1] * t + 0.5 * a[1] * t * t,
        p0[2] + v[2] * t + 0.5 * a[2] * t * t,
    )


def make_turn(omega=0.4, radius=120.0):
    return lambda t: (radius * math.cos(omega * t), radius * math.sin(omega * t), 8.0 + 0.5 * t)


def main():
    cases = [
        ("constant velocity", make_cv()),
        ("constant acceleration", make_ca()),
        ("coordinated turn", make_turn()),
    ]
    print(f"lead tau=1.0s, meas noise=0.5m, 20 Hz")
    print(f"{'track':<24} {'KF RMS (m)':>12} {'CV RMS (m)':>12} {'KF/CV':>8}")
    for name, track in cases:
        kf, cv = evaluate(track)
        print(f"{name:<24} {kf:>12.3f} {cv:>12.3f} {kf / cv:>8.2f}")


if __name__ == "__main__":
    main()
