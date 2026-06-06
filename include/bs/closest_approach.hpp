#pragma once

#include "target.hpp"

// ================================================================
// Golden-section on [0, dt]
// ================================================================
template <class EvalF>
inline double golden_section_min(EvalF&& f, double dt, int maxIter, double tolAbs, double tolRel)
{
    // Golden ratio constants
    constexpr double INVPHI = 0.6180339887498948482;
    constexpr double INVPHI2 = 0.3819660112501051518;

    double a = 0.0;
    double b = dt;

    double x1 = a + INVPHI2 * (b - a);
    double x2 = a + INVPHI * (b - a);

    double f1 = f(x1);
    double f2 = f(x2);

    const double tolTau = std::max(tolAbs, tolRel * dt);

    for (int k = 0; k < maxIter; ++k)
    {
        if ((b - a) <= tolTau)
        {
            break;
        }

        if (f2 < f1)
        {
            a = x1;
            x1 = x2;
            f1 = f2;

            x2 = a + INVPHI * (b - a);
            f2 = f(x2);
        }
        else
        {
            b = x2;
            x2 = x1;
            f2 = f1;

            x1 = a + INVPHI2 * (b - a);
            f1 = f(x1);
        }
    }

    double bestX = (f1 <= f2 ? x1 : x2);
    double bestF = (f1 <= f2 ? f1 : f2);

    const double f0 = f(0.0);
    if (f0 < bestF) { bestX = 0.0; bestF = f0; }

    const double fdt = f(dt);
    if (fdt < bestF) { bestX = dt; }

    return bestX;
}

// ================================================================
// Find closest approach (projectile starts at origin)
// ================================================================
inline void find_closest_approach_acc(
    const Vec3& projVel0,
    const Vec3& relPos0,
    const Vec3& relVel,
    const Vec3& relAcc,
    double kDrag,
    const BallisticParams& P,
    Vec3& relMissAtStar,
    double& tStar);

inline void find_closest_approach(
    const Vec3& projVel0,
    const Vec3& relPos0,
    const Vec3& relVel,
    double kDrag,
    const BallisticParams& P,
    Vec3& relMissAtStar,
    double& tStar)
{
    const Vec3 zeroAcc = { 0.0, 0.0, 0.0 };
    find_closest_approach_acc(projVel0, relPos0, relVel, zeroAcc, kDrag, P, relMissAtStar, tStar);
}

inline void find_closest_approach_acc(
    const Vec3& projVel0,
    const Vec3& relPos0,
    const Vec3& relVel,
    const Vec3& relAcc,
    double kDrag,
    const BallisticParams& P,
    Vec3& relMissAtStar,
    double& tStar)
{
    State prev{};
    prev.r = { 0.0, 0.0, 0.0 };
    prev.v = projVel0;

    double tPrev = 0.0;

    const Vec3 rel0 = rel_vec_acc(prev.r, relPos0, relVel, relAcc, 0.0);
    const Vec3 relVel0 = prev.v - relVel;
    double qPrev = dot(rel0, relVel0);

    bool isLow = (P.arcMode == ArcMode::Low);

    while (tPrev < P.tMax)
    {
        State curr = prev;
        rk4_step(curr, P.dt, kDrag, P.g, P.wind);
        const double tCurr = tPrev + P.dt;

        const Vec3 relCurr = rel_vec_acc(curr.r, relPos0, relVel, relAcc, tCurr);

        const Vec3 relVelCurr = curr.v - (relVel + tCurr * relAcc);
        const double qCurr = dot(relCurr, relVelCurr);

        bool allowCheck = isLow || curr.v.z <= 0.0;

        if (allowCheck && (qPrev <= 0.0) && (qCurr >= 0.0))
        {
            auto eval_rel2 = [&](double tau) -> double
            {
                State st = prev;
                if (tau > 1e-15)
                {
                    rk4_step(st, tau, kDrag, P.g, P.wind);
                }

                const double tOut = tPrev + tau;
                const Vec3 rel = rel_vec_acc(st.r, relPos0, relVel, relAcc, tOut);
                return dot(rel, rel);
            };

            const double tauStar = golden_section_min(eval_rel2, P.dt, P.gsMaxIter, P.gsTolAbs, P.gsTolRel);
            tStar = tPrev + tauStar;

            State stStar = prev;
            if (tauStar > 1e-15)
            {
                rk4_step(stStar, tauStar, kDrag, P.g, P.wind);
            }

            if (isLow || stStar.v.z <= 0.0)
            {
                relMissAtStar = rel_vec_acc(stStar.r, relPos0, relVel, relAcc, tStar);
                return;
            }
        }

        prev = curr;
        tPrev = tCurr;

        qPrev = qCurr;
    }

    tStar = tPrev;
    relMissAtStar = rel_vec_acc(prev.r, relPos0, relVel, relAcc, tStar);
}
