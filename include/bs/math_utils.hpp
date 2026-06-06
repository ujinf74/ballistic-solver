#pragma once

#include "vec3.hpp"

// ================================================================
// Utils
// ================================================================
inline double wrap_pi(double a)
{
    a = std::fmod(a + M_PI, 2.0 * M_PI);
    if (a < 0.0)
    {
        a += 2.0 * M_PI;
    }
    return a - M_PI;
}

inline Vec3 angles_to_dir(double theta, double phi)
{
    const double ct = std::cos(theta);
    return { ct * std::cos(phi), ct * std::sin(phi), std::sin(theta) };
}

inline void vec_to_angles(const Vec3& r, double& theta, double& phi)
{
    const double n = norm(r);
    if (n < 1e-12)
    {
        theta = 0.0;
        phi = 0.0;
        return;
    }

    const Vec3 d = (1.0 / n) * r;
    theta = std::asin(std::clamp(d.z, -1.0, 1.0));
    phi = std::atan2(d.y, d.x);
}

enum class ArcMode
{
    Low,
    High
};

enum class ParamPreset
{
    Fast = 0,
    Balanced = 1,
    Precise = 2
};
