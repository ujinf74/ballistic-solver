#pragma once

#include "integrator.hpp"

// ================================================================
// Relative geometry helpers
// ================================================================
inline Vec3 target_pos(const Vec3& relPos0, const Vec3& relVel, double t)
{
    return relPos0 + t * relVel;
}

inline Vec3 target_pos_acc(const Vec3& relPos0, const Vec3& relVel, const Vec3& relAcc, double t)
{
    return relPos0 + t * relVel + (0.5 * t * t) * relAcc;
}

inline Vec3 rel_vec(const Vec3& projPos, const Vec3& relPos0, const Vec3& relVel, double t)
{
    return projPos - target_pos(relPos0, relVel, t);
}

inline Vec3 rel_vec_acc(const Vec3& projPos, const Vec3& relPos0, const Vec3& relVel, const Vec3& relAcc, double t)
{
    return projPos - target_pos_acc(relPos0, relVel, relAcc, t);
}
