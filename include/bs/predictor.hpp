#pragma once

#include <functional>

#include "solve.hpp"
#include "solve_coord_lead.hpp"

// ================================================================
// Target predictor seam
// ================================================================
// A caller-supplied target trajectory in the projectile-relative frame:
// `pos(t)` returns the predicted relative position at time t >= 0. Velocity and
// acceleration are obtained by finite differences, so a position-only predictor
// (e.g. a Kalman/IMM track) is sufficient and relVel/relAcc need not be passed.
struct TargetPredictor
{
    std::function<Vec3(double)> pos;
    double fdStep = 1e-3;

    Vec3 position(double t) const
    {
        return pos(t);
    }

    Vec3 velocity(double t) const
    {
        const double h = fdStep;
        if (t - h < 0.0)
        {
            return (1.0 / h) * (pos(t + h) - pos(t));
        }
        return (1.0 / (2.0 * h)) * (pos(t + h) - pos(t - h));
    }

    Vec3 acceleration(double t) const
    {
        const double h = fdStep;
        if (t - h < 0.0)
        {
            return (1.0 / (h * h)) * (pos(t + 2.0 * h) - 2.0 * pos(t + h) + pos(t));
        }
        return (1.0 / (h * h)) * (pos(t + h) - 2.0 * pos(t) + pos(t - h));
    }
};

// ================================================================
// Solve against an arbitrary predictor (additive; core solver unchanged)
// ================================================================
// Reuses solve_launch_angles by fitting, at the current intercept-time estimate
// tStar, a local constant-acceleration model matching the predictor's position,
// velocity and acceleration there, then iterating tStar to a fixed point. This
// is exact for constant-velocity / constant-acceleration predictors and is
// second-order accurate at the intercept for smoothly curving trajectories.
inline SolverResult solve_launch_angles_predicted(
    const TargetPredictor& target,
    double v0,
    double kDrag,
    const BallisticParams& P = BallisticParams{},
    int maxOuter = 8,
    double tStarTol = 1e-3)
{
    SolverResult out{};

    double tStar = norm(target.position(0.0)) / std::max(v0, 1e-6);
    if (!std::isfinite(tStar) || tStar <= 0.0)
    {
        tStar = std::max(P.dt, 0.1 * P.tMax);
    }

    for (int k = 0; k < std::max(1, maxOuter); ++k)
    {
        const double tf = std::clamp(tStar, 0.0, P.tMax);

        const Vec3 p = target.position(tf);
        const Vec3 v = target.velocity(tf);
        const Vec3 a = target.acceleration(tf);

        // Local CA model expressed from t = 0 so that, at t = tf, it matches the
        // predictor's position/velocity/acceleration:
        const Vec3 relAcc = a;
        const Vec3 relVel = v - tf * a;
        const Vec3 relPos0 = p - tf * v + (0.5 * tf * tf) * a;

        out = solve_launch_angles_coord_lead(relPos0, relVel, v0, kDrag, P, relAcc);

        if (!std::isfinite(out.tStar))
        {
            break;
        }
        const double moved = std::fabs(out.tStar - tStar);
        tStar = out.tStar;
        if (moved <= tStarTol)
        {
            break;
        }
    }

    return out;
}
