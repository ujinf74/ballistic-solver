#pragma once

#include "params.hpp"

// ================================================================
// RK4 + quadratic drag
// ================================================================
struct State
{
    Vec3 r;
    Vec3 v;
};

inline void deriv(const State& s, double kDrag, double g, const Vec3& wind, State& ds)
{
    ds.r = s.v;

    Vec3 a = { 0.0, 0.0, -g };

    const Vec3 vRelAir = s.v - wind;

    const double speed = norm(vRelAir);
    if (speed > 1e-12 && kDrag != 0.0)
    {
        a = a - (kDrag * speed) * vRelAir;
    }

    ds.v = a;
}

inline void rk4_step(State& s, double h, double kDrag, double g, const Vec3& wind)
{
    State k1, k2, k3, k4, tmp;

    deriv(s, kDrag, g, wind, k1);

    tmp = s;
    tmp.r += (0.5 * h) * k1.r;
    tmp.v += (0.5 * h) * k1.v;
    deriv(tmp, kDrag, g, wind, k2);

    tmp = s;
    tmp.r += (0.5 * h) * k2.r;
    tmp.v += (0.5 * h) * k2.v;
    deriv(tmp, kDrag, g, wind, k3);

    tmp = s;
    tmp.r += h * k3.r;
    tmp.v += h * k3.v;
    deriv(tmp, kDrag, g, wind, k4);

    s.r += (h / 6.0) * (k1.r + 2.0 * k2.r + 2.0 * k3.r + k4.r);
    s.v += (h / 6.0) * (k1.v + 2.0 * k2.v + 2.0 * k3.v + k4.v);
}
