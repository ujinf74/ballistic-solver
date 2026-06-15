#pragma once

#include "solve.hpp"

// ================================================================
// Coordinate-residual launch solver (production).
//
// Solves launch angles by Gauss-Newton on the 3D closest-approach miss vector,
// with the initial inverse-Jacobian seeded ANALYTICALLY from the vacuum-arc map
// (no trajectory integration for the Jacobian) and refined by Broyden. No
// Levenberg-Marquardt damping, no line search. If the warm start does not reach
// tolerance, a small arc-appropriate theta-grid multistart is tried.
//
// Versus the auxiliary-residual solver (solve_launch_angles): equal robustness
// with the multistart fallback, ~15-55% faster, and simpler. The vacuum map is
// used as the seed and as the Jacobian preconditioner, not as the residual.
// ================================================================

// vacuum analytic inverse Jacobian M (2x3) = d(theta,phi)/d(aim), via finite
// differences of the closed-form vacuum_arc map -- NO trajectory integration.
inline bool vacuum_inverse_jacobian(const Vec3& aim, double v0, const BallisticParams& P, double M[2][3])
{
    double t0, p0;
    if (!vacuum_arc_angles_to_point(aim, v0, P.arcMode, P.g, t0, p0))
    {
        return false;
    }
    const double eps = std::clamp(1e-3 * (1.0 + norm(aim)), 1e-4, 1.0);
    for (int j = 0; j < 3; ++j)
    {
        Vec3 ap = aim, am = aim;
        if (j == 0) { ap.x += eps; am.x -= eps; }
        else if (j == 1) { ap.y += eps; am.y -= eps; }
        else { ap.z += eps; am.z -= eps; }
        double tp, pp, tm, pm;
        if (!vacuum_arc_angles_to_point(ap, v0, P.arcMode, P.g, tp, pp) ||
            !vacuum_arc_angles_to_point(am, v0, P.arcMode, P.g, tm, pm))
        {
            return false;
        }
        M[0][j] = (tp - tm) / (2.0 * eps);
        M[1][j] = wrap_pi(pp - pm) / (2.0 * eps);
    }
    return std::isfinite(M[0][0]) && std::isfinite(M[1][2]);
}

// One coordinate-residual solve from a fixed (theta, phi) seed.
inline SolverResult coord_vac_once(
    const Vec3& relPos0, const Vec3& relVel, const Vec3& relAcc,
    double v0, double kDrag, const BallisticParams& P,
    double seedTheta, double seedPhi)
{
    auto residual = [&](double th, double ph, Vec3& r, double& miss, double& tStar)
    {
        find_closest_approach_acc(v0 * angles_to_dir(th, ph), relPos0, relVel, relAcc, kDrag, P, r, tStar);
        miss = norm(r);
    };

    double th = std::clamp(seedTheta, P.thetaMin, P.thetaMax);
    double ph = wrap_pi(seedPhi);
    Vec3 r{};
    double miss = std::numeric_limits<double>::quiet_NaN();
    double tStar = std::numeric_limits<double>::quiet_NaN();
    residual(th, ph, r, miss, tStar);

    double bTh = th, bPh = ph, bMiss = miss, bTime = tStar;
    Vec3 bRel = r;

    const Vec3 aim = target_pos_acc(relPos0, relVel, relAcc, tStar);
    double M[2][3];
    if (!vacuum_inverse_jacobian(aim, v0, P, M))
    {
        SolverResult out{};
        out.theta = bTh; out.phi = bPh; out.miss = bMiss; out.relMissAtStar = bRel; out.tStar = bTime;
        out.success = std::isfinite(bMiss) && bMiss <= P.tolMiss;
        out.report.status = out.success ? SolveStatus::Ok : SolveStatus::JacobianFailed;
        return out;
    }

    int it = 0;
    for (; it < P.maxIter; ++it)
    {
        if (miss <= P.tolMiss)
        {
            break;
        }
        const double dth = -(M[0][0] * r.x + M[0][1] * r.y + M[0][2] * r.z);
        const double dph = -(M[1][0] * r.x + M[1][1] * r.y + M[1][2] * r.z);
        const double nth = std::clamp(th + dth, P.thetaMin, P.thetaMax);
        const double nph = wrap_pi(ph + dph);
        Vec3 rn{};
        double mn, tn;
        residual(nth, nph, rn, mn, tn);

        // bad-Broyden inverse update: M += ((s - M y) y^T) / (y^T y)
        const double s0 = nth - th, s1 = wrap_pi(nph - ph);
        const Vec3 y = rn - r;
        const double yy = dot(y, y);
        if (yy > 1e-18)
        {
            const double My0 = M[0][0] * y.x + M[0][1] * y.y + M[0][2] * y.z;
            const double My1 = M[1][0] * y.x + M[1][1] * y.y + M[1][2] * y.z;
            const double c0 = (s0 - My0) / yy, c1 = (s1 - My1) / yy;
            M[0][0] += c0 * y.x; M[0][1] += c0 * y.y; M[0][2] += c0 * y.z;
            M[1][0] += c1 * y.x; M[1][1] += c1 * y.y; M[1][2] += c1 * y.z;
        }

        th = nth; ph = nph; r = rn; miss = mn; tStar = tn;
        if (std::isfinite(miss) && miss < bMiss)
        {
            bTh = th; bPh = ph; bMiss = miss; bRel = r; bTime = tStar;
        }
    }

    SolverResult out{};
    out.theta = bTh; out.phi = bPh; out.miss = bMiss; out.relMissAtStar = bRel; out.tStar = bTime;
    out.success = std::isfinite(bMiss) && bMiss <= P.tolMiss;
    out.report.iterations = it;
    out.report.lastMiss = bMiss;
    out.report.status = out.success ? SolveStatus::Ok : SolveStatus::MaxIterReached;
    out.report.message = out.success ? "Ok" : "MaxIterReached: coordinate-residual solver did not reach tolMiss.";
    return out;
}

// Public entry: warm start from the vacuum lead, then arc-grid multistart fallback.
inline SolverResult solve_launch_angles_coord_lead(
    const Vec3& relPos0,
    const Vec3& relVel,
    double v0,
    double kDrag,
    const BallisticParams& P = BallisticParams{},
    const Vec3& relAcc = Vec3{ 0.0, 0.0, 0.0 })
{
    SolverResult out{};
    if (!std::isfinite(v0) || v0 <= 0.0 || !std::isfinite(kDrag) ||
        !std::isfinite(P.g) || P.g <= 0.0 || !std::isfinite(P.dt) || P.dt <= 0.0 ||
        !std::isfinite(P.tMax) || P.tMax <= 0.0 || P.maxIter <= 0 ||
        P.thetaMin >= P.thetaMax)
    {
        out.report.status = SolveStatus::InvalidInput;
        out.report.message = "InvalidInput: v0/g/dt/tMax/maxIter/theta range check failed.";
        return out;
    }

    double seedTh, seedPh;
    initial_guess_vacuum_lead_acc(relPos0, relVel, relAcc, v0, P.arcMode, P.g, seedTh, seedPh, P.tMax);
    seedTh = std::clamp(seedTh, P.thetaMin, P.thetaMax);
    seedPh = wrap_pi(seedPh);

    SolverResult best = coord_vac_once(relPos0, relVel, relAcc, v0, kDrag, P, seedTh, seedPh);
    if (best.success)
    {
        return best;
    }

    // Multistart fallback: arc-appropriate theta grid, azimuth from the vacuum seed.
    static const double kHigh[] = { 50.0, 54.0, 58.0, 62.0, 66.0, 70.0, 74.0, 78.0 };
    static const double kLow[] = { 6.0, 12.0, 18.0, 24.0, 30.0, 36.0 };
    const double* grid = (P.arcMode == ArcMode::High) ? kHigh : kLow;
    const int nGrid = (P.arcMode == ArcMode::High) ? 8 : 6;
    for (int i = 0; i < nGrid; ++i)
    {
        const double th = std::clamp(grid[i] * (M_PI / 180.0), P.thetaMin, P.thetaMax);
        SolverResult r = coord_vac_once(relPos0, relVel, relAcc, v0, kDrag, P, th, seedPh);
        if (std::isfinite(r.miss) && r.miss < best.miss)
        {
            best = r;
        }
        if (r.success)
        {
            return r;
        }
    }
    return best;
}
