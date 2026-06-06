#pragma once

#include "closest_approach.hpp"

// ================================================================
// Vacuum_arc_angles_to_point
// ================================================================
inline bool vacuum_arc_angles_to_point(const Vec3& R, double v0, ArcMode mode, double g, double& theta, double& phi)
{
    const double Rxy = std::sqrt(R.x * R.x + R.y * R.y);
    phi = std::atan2(R.y, R.x);

    if (Rxy < 1e-12)
    {
        theta = (R.z >= 0.0 ? 1.0 : -1.0) * (M_PI / 2.0);
        return true;
    }

    const double v2 = v0 * v0;
    const double disc = v2 * v2 - g * (g * Rxy * Rxy + 2.0 * R.z * v2);

    if (disc < 0.0)
    {
        theta = std::atan2(R.z, Rxy);
        return false;
    }

    const double s = std::sqrt(disc);
    const double num = (mode == ArcMode::High) ? (v2 + s) : (v2 - s);
    const double tanTheta = num / (g * Rxy);

    theta = std::atan(tanTheta);

    if (!std::isfinite(theta))
    {
        theta = std::atan2(R.z, Rxy);
        return false;
    }

    return true;
}

// ================================================================
// Initial guess (vacuum lead)
// ================================================================
inline void initial_guess_vacuum_lead_acc(
    const Vec3& relPos0,
    const Vec3& relVel,
    const Vec3& relAcc,
    double v0,
    ArcMode mode,
    double g,
    double& theta0,
    double& phi0,
    double tMax = 40.0);

inline void initial_guess_vacuum_lead(
    const Vec3& relPos0,
    const Vec3& relVel,
    double v0,
    ArcMode mode,
    double g,
    double& theta0,
    double& phi0)
{
    const Vec3 zeroAcc = { 0.0, 0.0, 0.0 };
    initial_guess_vacuum_lead_acc(relPos0, relVel, zeroAcc, v0, mode, g, theta0, phi0);
}

inline double vacuum_lead_linear_time_residual(
    const Vec3& relPos0,
    const Vec3& relVel,
    double v0,
    double g,
    double t)
{
    const Vec3 aim = relPos0 + t * relVel + Vec3{ 0.0, 0.0, 0.5 * g * t * t };
    return dot(aim, aim) - v0 * v0 * t * t;
}

inline double eval_quartic(double a4, double a3, double a2, double a1, double a0, double x)
{
    return (((a4 * x + a3) * x + a2) * x + a1) * x + a0;
}

inline double eval_quartic_deriv(double a4, double a3, double a2, double a1, double x)
{
    return ((4.0 * a4 * x + 3.0 * a3) * x + 2.0 * a2) * x + a1;
}

inline double polish_quartic_root(double a4, double a3, double a2, double a1, double a0, double x)
{
    for (int i = 0; i < 8; ++i)
    {
        const double f = eval_quartic(a4, a3, a2, a1, a0, x);
        const double df = eval_quartic_deriv(a4, a3, a2, a1, x);
        if (!std::isfinite(f) || !std::isfinite(df) || std::fabs(df) < 1e-14)
        {
            break;
        }
        const double dx = f / df;
        x -= dx;
        if (!std::isfinite(x) || std::fabs(dx) <= 1e-12 * (1.0 + std::fabs(x)))
        {
            break;
        }
    }
    return x;
}

inline bool add_quartic_real_root(double roots[4], int& count, double root)
{
    if (!std::isfinite(root))
    {
        return false;
    }

    for (int i = 0; i < count; ++i)
    {
        if (std::fabs(roots[i] - root) <= 1e-7 * (1.0 + std::max(std::fabs(roots[i]), std::fabs(root))))
        {
            return false;
        }
    }

    if (count >= 4)
    {
        return false;
    }

    roots[count++] = root;
    return true;
}

inline int solve_quartic_real_roots(
    double a4,
    double a3,
    double a2,
    double a1,
    double a0,
    double roots[4])
{
    if (!std::isfinite(a4) || !std::isfinite(a3) || !std::isfinite(a2) ||
        !std::isfinite(a1) || !std::isfinite(a0) || std::fabs(a4) < 1e-30)
    {
        return 0;
    }

    using Complex = std::complex<double>;

    const double A = a3 / a4;
    const double B = a2 / a4;
    const double C = a1 / a4;
    const double D = a0 / a4;

    const double p = B - 3.0 * A * A / 8.0;
    const double q = A * A * A / 8.0 - A * B / 2.0 + C;
    const double r = -3.0 * A * A * A * A / 256.0 + A * A * B / 16.0 - A * C / 4.0 + D;

    std::array<Complex, 4> candidates{};
    int candidateCount = 0;

    if (std::fabs(q) <= 1e-14 * (1.0 + std::fabs(p) + std::fabs(r)))
    {
        const Complex disc = Complex(p * p - 4.0 * r, 0.0);
        const Complex z0 = (-p + std::sqrt(disc)) * 0.5;
        const Complex z1 = (-p - std::sqrt(disc)) * 0.5;
        candidates[candidateCount++] = std::sqrt(z0) - A / 4.0;
        candidates[candidateCount++] = -std::sqrt(z0) - A / 4.0;
        candidates[candidateCount++] = std::sqrt(z1) - A / 4.0;
        candidates[candidateCount++] = -std::sqrt(z1) - A / 4.0;
    }
    else
    {
        const Complex Pc = -p * p / 12.0 - r;
        const Complex Qc = -p * p * p / 108.0 + p * r / 3.0 - q * q / 8.0;
        const Complex sqrtArg = Qc * Qc / 4.0 + Pc * Pc * Pc / 27.0;
        Complex U = std::pow(-Qc / 2.0 + std::sqrt(sqrtArg), 1.0 / 3.0);
        if (std::abs(U) <= 1e-14)
        {
            U = std::pow(-Qc / 2.0 - std::sqrt(sqrtArg), 1.0 / 3.0);
        }

        if (std::abs(U) > 1e-14)
        {
            const Complex y = -5.0 * p / 6.0 + U - Pc / (3.0 * U);
            const Complex W = std::sqrt(p + 2.0 * y);
            if (std::abs(W) > 1e-14)
            {
                const Complex common = -3.0 * p - 2.0 * y;
                const Complex term0 = std::sqrt(common - 2.0 * q / W);
                const Complex term1 = std::sqrt(common + 2.0 * q / W);
                candidates[candidateCount++] = 0.5 * ( W + term0) - A / 4.0;
                candidates[candidateCount++] = 0.5 * ( W - term0) - A / 4.0;
                candidates[candidateCount++] = 0.5 * (-W + term1) - A / 4.0;
                candidates[candidateCount++] = 0.5 * (-W - term1) - A / 4.0;
            }
        }
    }

    int count = 0;
    for (int i = 0; i < candidateCount; ++i)
    {
        const double real = candidates[i].real();
        const double imag = candidates[i].imag();
        if (std::fabs(imag) <= 1e-6 * (1.0 + std::fabs(real)))
        {
            const double polished = polish_quartic_root(a4, a3, a2, a1, a0, real);
            add_quartic_real_root(roots, count, polished);
        }
    }

    return count;
}

inline bool vacuum_lead_quartic_linear_target(
    const Vec3& relPos0,
    const Vec3& relVel,
    double v0,
    ArcMode mode,
    double g,
    double tMax,
    double& theta,
    double& phi)
{
    if (!std::isfinite(v0) || v0 <= 0.0 || !std::isfinite(g) || g <= 0.0 ||
        !std::isfinite(tMax) || tMax <= 1e-6)
    {
        return false;
    }

    bool have = false;
    double bestTheta = 0.0;
    double bestPhi = 0.0;

    const Vec3 q = { 0.0, 0.0, 0.5 * g };
    const double a4 = dot(q, q);
    const double a3 = 2.0 * dot(relVel, q);
    const double a2 = dot(relVel, relVel) + 2.0 * dot(relPos0, q) - v0 * v0;
    const double a1 = 2.0 * dot(relPos0, relVel);
    const double a0 = dot(relPos0, relPos0);

    double roots[4];
    const int rootCount = solve_quartic_real_roots(a4, a3, a2, a1, a0, roots);
    for (int i = 0; i < rootCount; ++i)
    {
        const double root = roots[i];
        if (root <= 1e-8 || root > tMax * (1.0 + 1e-9))
        {
            continue;
        }

        const Vec3 aim = relPos0 + root * relVel + Vec3{ 0.0, 0.0, 0.5 * g * root * root };
        double th;
        double ph;
        vec_to_angles(aim, th, ph);

        if (std::isfinite(th) && std::isfinite(ph) &&
            (!have || (mode == ArcMode::High ? th > bestTheta : th < bestTheta)))
        {
            have = true;
            bestTheta = th;
            bestPhi = ph;
        }
    }

    if (!have)
    {
        return false;
    }

    theta = bestTheta;
    phi = bestPhi;
    return true;
}

inline void initial_guess_vacuum_lead_acc(
    const Vec3& relPos0,
    const Vec3& relVel,
    const Vec3& relAcc,
    double v0,
    ArcMode mode,
    double g,
    double& theta0,
    double& phi0,
    double tMax)
{
    if (norm(relAcc) <= 1e-12 &&
        vacuum_lead_quartic_linear_target(relPos0, relVel, v0, mode, g, tMax, theta0, phi0))
    {
        return;
    }

    const double t0 = norm(relPos0) / std::max(v0, 1e-6);
    const Vec3 R = target_pos_acc(relPos0, relVel, relAcc, t0);

    const double Rxy = std::sqrt(R.x * R.x + R.y * R.y);
    phi0 = std::atan2(R.y, R.x);

    if (Rxy < 1e-9)
    {
        theta0 = (R.z >= 0.0 ? 1.0 : -1.0) * (M_PI / 2.0);
        return;
    }

    const double v2 = v0 * v0;
    const double disc = v2 * v2 - g * (g * Rxy * Rxy + 2.0 * R.z * v2);

    if (disc < 0.0)
    {
        theta0 = std::atan2(R.z, Rxy);
        return;
    }

    const double s = std::sqrt(disc);
    const double num = (mode == ArcMode::High) ? (v2 + s) : (v2 - s);
    const double tanTheta = num / (g * Rxy);

    theta0 = std::atan(tanTheta);

    if (!std::isfinite(theta0))
    {
        const double tanHigh = (v2 + std::sqrt(disc)) / (g * Rxy);
        theta0 = std::atan(tanHigh);
    }
}
