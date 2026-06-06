#pragma once

#include "vacuum_lead.hpp"

// ================================================================
// Residual
// ================================================================
inline bool compute_angle_residual_acc(
    double theta,
    double phi,
    const Vec3& relPos0,
    const Vec3& relVel,
    const Vec3& relAcc,
    double v0,
    double kDrag,
    const BallisticParams& P,
    double F[2],
    double& miss,
    Vec3& relMissAtStar_out,
    double& tStar_out);

inline bool compute_angle_residual(
    double theta,
    double phi,
    const Vec3& relPos0,
    const Vec3& relVel,
    double v0,
    double kDrag,
    const BallisticParams& P,
    double F[2],
    double& miss,
    Vec3& relMissAtStar_out,
    double& tStar_out)
{
    const Vec3 zeroAcc = { 0.0, 0.0, 0.0 };
    return compute_angle_residual_acc(theta, phi, relPos0, relVel, zeroAcc, v0, kDrag, P, F, miss, relMissAtStar_out, tStar_out);
}

inline bool compute_auxiliary_delta(
    const Vec3& aim,
    const Vec3& relMissAtStar,
    double beta,
    double v0,
    const BallisticParams& P,
    double& dtheta,
    double& dphi)
{
    const Vec3 aimCorr = aim - beta * relMissAtStar;

    double th0, ph0, th1, ph1;
    const bool ok0 = vacuum_arc_angles_to_point(aim, v0, P.arcMode, P.g, th0, ph0);
    const bool ok1 = vacuum_arc_angles_to_point(aimCorr, v0, P.arcMode, P.g, th1, ph1);

    if (!ok0 || !ok1)
    {
        vec_to_angles(aim, th0, ph0);
        vec_to_angles(aimCorr, th1, ph1);
    }

    dtheta = th1 - th0;
    dphi = wrap_pi(ph1 - ph0);

    return std::isfinite(dtheta) && std::isfinite(dphi);
}

inline bool compute_angle_residual_acc(
    double theta,
    double phi,
    const Vec3& relPos0,
    const Vec3& relVel,
    const Vec3& relAcc,
    double v0,
    double kDrag,
    const BallisticParams& P,
    double F[2],
    double& miss,
    Vec3& relMissAtStar_out,
    double& tStar_out)
{
    const Vec3 dir = angles_to_dir(theta, phi);
    const Vec3 projVel0 = v0 * dir;

    Vec3 relMissAtStar;
    double tStar;
    find_closest_approach_acc(projVel0, relPos0, relVel, relAcc, kDrag, P, relMissAtStar, tStar);

    miss = norm(relMissAtStar);
    relMissAtStar_out = relMissAtStar;
    tStar_out = tStar;

    const Vec3 aim = target_pos_acc(relPos0, relVel, relAcc, tStar);
    double dtheta, dphi;
    if (!compute_auxiliary_delta(aim, relMissAtStar, P.beta, v0, P, dtheta, dphi))
    {
        return false;
    }

    F[0] = dtheta;
    F[1] = dphi;

    return std::isfinite(miss);
}
