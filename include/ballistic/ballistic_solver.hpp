#pragma once

#include <cmath>
#include <algorithm>
#include <limits>
#include <string>

#ifndef M_PI
#define M_PI 3.141592653589793238462643383279502884
#endif

// ================================================================
// Vec3 (minimal)
// ================================================================
struct Vec3
{
    double x, y, z;
};

inline Vec3 operator+(const Vec3& a, const Vec3& b)
{
    return { a.x + b.x, a.y + b.y, a.z + b.z };
}

inline Vec3 operator-(const Vec3& a, const Vec3& b)
{
    return { a.x - b.x, a.y - b.y, a.z - b.z };
}

inline Vec3 operator*(double s, const Vec3& v)
{
    return { s * v.x, s * v.y, s * v.z };
}

inline Vec3& operator+=(Vec3& a, const Vec3& b)
{
    a.x += b.x;
    a.y += b.y;
    a.z += b.z;
    return a;
}

inline double dot(const Vec3& a, const Vec3& b)
{
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

inline double norm(const Vec3& v)
{
    return std::sqrt(dot(v, v));
}

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

// ================================================================
// // Parameters (configurable)
// ================================================================
struct BallisticParams
{
    // ----------------------------
    // Physical constants
    // ----------------------------
    double g = 9.80665;        // Gravity acceleration (m/s^2)

    // ----------------------------
    // Integration/search parameters
    // ----------------------------
    double dt = 0.01;          // RK4 basic step
    double tMax = 20.0;        // Maximum simulation time

    // ----------------------------
    // Targeting/acceptance criteria
    // ----------------------------
    double tolMiss = 1e-2;     // Success tolerance (miss distance)
    double beta = 1.0;         // Target offset correction factor

    // ----------------------------
    // Solver iteration control
    // ----------------------------
    int maxIter = 20;

    // LM(lambda) control
    double lambdaInit = 1e-6;
    double lambdaMin = 1e-12;
    double lambdaMax = 1e+6;
    double lambdaUpMul = 10.0;
    double lambdaDownMul = 0.3;
    int lambdaTries = 6;

    // Line search
    int lineSearchTries = 10;
    double alphaMin = 1e-4;

    // Angle bounds
    double thetaMin = -M_PI / 2.0;
    double thetaMax =  M_PI / 2.0;

    // Broyden
    double broydenMinDenom = 1e-12;

    // Acceptance tolerance relaxation
    double missEps = 1e-12;

    // FD step size
    double fdScale = 1e-4;
    double fdMin = 1e-6;
    double fdMax = 1e-3;

    // Golden-section
    int gsMaxIter = 20;
    double gsTolAbs = 1e-8;
    double gsTolRel = 1e-8;
};

// ================================================================
// Error reporting
// ================================================================
enum class SolveStatus
{
    Ok = 0,

    // Input / evaluation stage
    InvalidInput,
    InitialResidualFailed,
    JacobianFailed,

    // During iteration
    LMStepSingular,
    ResidualFailedDuringSearch,
    LineSearchRejected,
    LambdaTriesExhausted,
    MaxIterReached
};

struct SolveReport
{
    SolveStatus status = SolveStatus::Ok;
    std::string message;

    // For debugging / diagnostics
    int iterations = 0;
    int acceptedSteps = 0;

    double lastLambda = std::numeric_limits<double>::quiet_NaN();
    double lastAlpha = std::numeric_limits<double>::quiet_NaN();

    double lastTheta = std::numeric_limits<double>::quiet_NaN();
    double lastPhi = std::numeric_limits<double>::quiet_NaN();

    double lastMiss = std::numeric_limits<double>::quiet_NaN();
    double lastF0 = std::numeric_limits<double>::quiet_NaN();
    double lastF1 = std::numeric_limits<double>::quiet_NaN();
};

// ================================================================
// RK4 + quadratic drag
// ================================================================
struct State
{
    Vec3 r;
    Vec3 v;
};

inline void deriv(const State& s, double kDrag, double g, State& ds)
{
    ds.r = s.v;

    Vec3 a = { 0.0, 0.0, -g };

    const double speed = norm(s.v);
    if (speed > 1e-12 && kDrag != 0.0)
    {
        a = a - (kDrag * speed) * s.v;
    }

    ds.v = a;
}

inline void rk4_step(State& s, double h, double kDrag, double g)
{
    State k1, k2, k3, k4, tmp;

    deriv(s, kDrag, g, k1);

    tmp = s;
    tmp.r += (0.5 * h) * k1.r;
    tmp.v += (0.5 * h) * k1.v;
    deriv(tmp, kDrag, g, k2);

    tmp = s;
    tmp.r += (0.5 * h) * k2.r;
    tmp.v += (0.5 * h) * k2.v;
    deriv(tmp, kDrag, g, k3);

    tmp = s;
    tmp.r += h * k3.r;
    tmp.v += h * k3.v;
    deriv(tmp, kDrag, g, k4);

    s.r += (h / 6.0) * (k1.r + 2.0 * k2.r + 2.0 * k3.r + k4.r);
    s.v += (h / 6.0) * (k1.v + 2.0 * k2.v + 2.0 * k3.v + k4.v);
}

// ================================================================
// Relative geometry helpers
// ================================================================
inline Vec3 target_pos(const Vec3& relPos0, const Vec3& relVel, double t)
{
    return relPos0 + t * relVel;
}

inline Vec3 rel_vec(const Vec3& projPos, const Vec3& relPos0, const Vec3& relVel, double t)
{
    return projPos - target_pos(relPos0, relVel, t);
}

inline bool has_turn(bool hadDecreased, double d2Prev, double d2Curr, double qPrev, double qCurr)
{
    const bool turnByDistance = (hadDecreased && (d2Curr > d2Prev));
    const bool turnByQSign = (hadDecreased && (qPrev <= 0.0) && (qCurr >= 0.0));
    return turnByDistance || turnByQSign;
}

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

    return std::clamp((f1 <= f2 ? x1 : x2), 0.0, dt);
}

// ================================================================
// Closest approach (projectile starts at origin)
// ================================================================
inline void closest_approach(
    const Vec3& projVel0,
    const Vec3& relPos0,
    const Vec3& relVel,
    double kDrag,
    const BallisticParams& P,
    Vec3& relMissAtStar,
    double& tStar)
{
    State prev{};
    prev.r = { 0.0, 0.0, 0.0 };
    prev.v = projVel0;

    double tPrev = 0.0;

    Vec3 relPrev = rel_vec(prev.r, relPos0, relVel, tPrev);
    double d2Prev = dot(relPrev, relPrev);

    Vec3 relVelPrev = prev.v - relVel;
    double qPrev = dot(relPrev, relVelPrev);

    bool hadDecreased = false;

    while (tPrev < P.tMax)
    {
        State curr = prev;
        rk4_step(curr, P.dt, kDrag, P.g);
        const double tCurr = tPrev + P.dt;

        const Vec3 relCurr = rel_vec(curr.r, relPos0, relVel, tCurr);
        const double d2Curr = dot(relCurr, relCurr);

        if (d2Curr < d2Prev)
        {
            hadDecreased = true;
        }

        const Vec3 relVelCurr = curr.v - relVel;
        const double qCurr = dot(relCurr, relVelCurr);

        if (has_turn(hadDecreased, d2Prev, d2Curr, qPrev, qCurr))
        {
            auto eval_rel2 = [&](double tau) -> double
            {
                State st = prev;
                if (tau > 1e-15)
                {
                    rk4_step(st, tau, kDrag, P.g);
                }

                const double tOut = tPrev + tau;
                const Vec3 rel = rel_vec(st.r, relPos0, relVel, tOut);
                return dot(rel, rel);
            };

            const double tauStar = golden_section_min(eval_rel2, P.dt, P.gsMaxIter, P.gsTolAbs, P.gsTolRel);
            tStar = tPrev + tauStar;

            State stStar = prev;
            if (tauStar > 1e-15)
            {
                rk4_step(stStar, tauStar, kDrag, P.g);
            }

            relMissAtStar = rel_vec(stStar.r, relPos0, relVel, tStar);
            return;
        }

        prev = curr;
        tPrev = tCurr;

        relPrev = relCurr;
        d2Prev = d2Curr;
        qPrev = qCurr;
    }

    tStar = tPrev;
    relMissAtStar = rel_vec(prev.r, relPos0, relVel, tStar);
}

// ================================================================
// Vacuum low-arc to point
// ================================================================
inline bool vacuum_low_arc_to_point(const Vec3& R, double v0, double g, double& theta, double& phi)
{
    const double Rxy = std::sqrt(R.x * R.x + R.y * R.y);
    phi = std::atan2(R.y, R.x);

    if (Rxy < 1e-12)
    {
        theta = (R.z >= 0.0 ? 1.0 : -1.0) * (M_PI / 4.0);
        return true;
    }

    const double v2 = v0 * v0;
    const double disc = v2 * v2 - g * (g * Rxy * Rxy + 2.0 * R.z * v2);

    if (disc < 0.0)
    {
        theta = std::atan2(R.z, Rxy);
        return false;
    }

    const double tanLow = (v2 - std::sqrt(disc)) / (g * Rxy);
    theta = std::atan(tanLow);

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
inline void initial_guess_vacuum_lead(
    const Vec3& relPos0,
    const Vec3& relVel,
    double v0,
    double g,
    double& theta0,
    double& phi0)
{
    const double t0 = norm(relPos0) / std::max(v0, 1e-6);
    const Vec3 R = relPos0 + t0 * relVel;

    const double Rxy = std::sqrt(R.x * R.x + R.y * R.y);
    phi0 = std::atan2(R.y, R.x);

    if (Rxy < 1e-9)
    {
        theta0 = (R.z >= 0.0 ? 1.0 : -1.0) * (M_PI / 4.0);
        return;
    }

    const double v2 = v0 * v0;
    const double disc = v2 * v2 - g * (g * Rxy * Rxy + 2.0 * R.z * v2);

    if (disc < 0.0)
    {
        theta0 = std::atan2(R.z, Rxy);
        return;
    }

    const double tanLow = (v2 - std::sqrt(disc)) / (g * Rxy);
    theta0 = std::atan(tanLow);

    if (!std::isfinite(theta0))
    {
        const double tanHigh = (v2 + std::sqrt(disc)) / (g * Rxy);
        theta0 = std::atan(tanHigh);
    }
}

// ================================================================
// Residual
// ================================================================
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
    const Vec3 dir = angles_to_dir(theta, phi);
    const Vec3 projVel0 = v0 * dir;

    Vec3 relMissAtStar;
    double tStar;
    closest_approach(projVel0, relPos0, relVel, kDrag, P, relMissAtStar, tStar);

    miss = norm(relMissAtStar);
    relMissAtStar_out = relMissAtStar;
    tStar_out = tStar;

    const Vec3 aim = target_pos(relPos0, relVel, tStar);
    const Vec3 aimCorr = aim - P.beta * relMissAtStar;

    double th0, ph0, th1, ph1;
    const bool ok0 = vacuum_low_arc_to_point(aim, v0, P.g, th0, ph0);
    const bool ok1 = vacuum_low_arc_to_point(aimCorr, v0, P.g, th1, ph1);

    if (!ok0 || !ok1)
    {
        vec_to_angles(aim, th0, ph0);
        vec_to_angles(aimCorr, th1, ph1);
    }

    F[0] = th1 - th0;
    F[1] = wrap_pi(ph1 - ph0);

    return std::isfinite(F[0]) && std::isfinite(F[1]) && std::isfinite(miss);
}

// ================================================================
// Solver bits
// ================================================================
inline bool solve_lm_step_2x2(const double J[2][2], const double F[2], double lambda, double& dtheta, double& dphi)
{
    const double A00 = J[0][0] * J[0][0] + J[1][0] * J[1][0] + lambda;
    const double A01 = J[0][0] * J[0][1] + J[1][0] * J[1][1];
    const double A11 = J[0][1] * J[0][1] + J[1][1] * J[1][1] + lambda;

    const double b0 = -(J[0][0] * F[0] + J[1][0] * F[1]);
    const double b1 = -(J[0][1] * F[0] + J[1][1] * F[1]);

    const double det = A00 * A11 - A01 * A01;
    if (std::fabs(det) < 1e-18)
    {
        return false;
    }

    dtheta = (b0 * A11 - b1 * A01) / det;
    dphi = (A00 * b1 - A01 * b0) / det;

    return std::isfinite(dtheta) && std::isfinite(dphi);
}

inline void broyden_update(double J[2][2], const double sdx[2], const double ydf[2], double minDenom)
{
    const double denom = sdx[0] * sdx[0] + sdx[1] * sdx[1];
    if (denom <= minDenom)
    {
        return;
    }

    const double Jsdx0 = J[0][0] * sdx[0] + J[0][1] * sdx[1];
    const double Jsdx1 = J[1][0] * sdx[0] + J[1][1] * sdx[1];

    const double u0 = (ydf[0] - Jsdx0) / denom;
    const double u1 = (ydf[1] - Jsdx1) / denom;

    J[0][0] += u0 * sdx[0];
    J[0][1] += u0 * sdx[1];
    J[1][0] += u1 * sdx[0];
    J[1][1] += u1 * sdx[1];
}

// ================================================================
// Result
// ================================================================
struct SolverResult
{
    bool success = false;
    double theta = std::numeric_limits<double>::quiet_NaN();
    double phi = std::numeric_limits<double>::quiet_NaN();
    double miss = std::numeric_limits<double>::quiet_NaN();
    Vec3 relMissAtStar = { std::numeric_limits<double>::quiet_NaN(),
                           std::numeric_limits<double>::quiet_NaN(),
                           std::numeric_limits<double>::quiet_NaN() };
    double tStar = std::numeric_limits<double>::quiet_NaN();

    SolveReport report;
};

// ================================================================
// Solve launch angles (Broyden + LM + line search)
// ================================================================
inline SolverResult solve_launch_angles(
    const Vec3& relPos0,
    const Vec3& relVel,
    double v0,
    double kDrag,
    const BallisticParams& P = BallisticParams{})
{
    SolverResult out{};

    // ----------------------------
    // Input validation
    // ----------------------------
    if (!std::isfinite(v0) || v0 <= 0.0 || !std::isfinite(kDrag) ||
        !std::isfinite(P.g) || P.g <= 0.0 || !std::isfinite(P.dt) || P.dt <= 0.0 ||
        !std::isfinite(P.tMax) || P.tMax <= 0.0 || P.maxIter <= 0 ||
        P.thetaMin >= P.thetaMax)
    {
        out.report.status = SolveStatus::InvalidInput;
        out.report.message = "InvalidInput: v0/g/dt/tMax/maxIter/theta range check failed.";
        return out;
    }

    auto pick_fd_step = [&](double x) -> double
    {
        double h = P.fdScale * (1.0 + std::fabs(x));
        return std::clamp(h, P.fdMin, P.fdMax);
    };

    auto jacobian_fd = [&](double th, double ph, const double Fbase[2], double Jout[2][2]) -> bool
    {
        double Fp[2], Fm[2];
        double mTmp;
        Vec3 relTmp;
        double tTmp;

        const double hth = pick_fd_step(th);

        const bool canMinus = (th - hth >= P.thetaMin);
        const bool canPlus = (th + hth <= P.thetaMax);

        if (canMinus && canPlus)
        {
            if (!compute_angle_residual(th + hth, ph, relPos0, relVel, v0, kDrag, P, Fp, mTmp, relTmp, tTmp))
            {
                return false;
            }
            if (!compute_angle_residual(th - hth, ph, relPos0, relVel, v0, kDrag, P, Fm, mTmp, relTmp, tTmp))
            {
                return false;
            }

            Jout[0][0] = (Fp[0] - Fm[0]) / (2.0 * hth);
            Jout[1][0] = wrap_pi(Fp[1] - Fm[1]) / (2.0 * hth);
        }
        else if (canPlus)
        {
            if (!compute_angle_residual(th + hth, ph, relPos0, relVel, v0, kDrag, P, Fp, mTmp, relTmp, tTmp))
            {
                return false;
            }

            Jout[0][0] = (Fp[0] - Fbase[0]) / hth;
            Jout[1][0] = wrap_pi(Fp[1] - Fbase[1]) / hth;
        }
        else if (canMinus)
        {
            if (!compute_angle_residual(th - hth, ph, relPos0, relVel, v0, kDrag, P, Fm, mTmp, relTmp, tTmp))
            {
                return false;
            }

            Jout[0][0] = (Fbase[0] - Fm[0]) / hth;
            Jout[1][0] = wrap_pi(Fbase[1] - Fm[1]) / hth;
        }
        else
        {
            return false;
        }

        const double hph = pick_fd_step(ph);

        const double php = wrap_pi(ph + hph);
        const double phm = wrap_pi(ph - hph);

        if (!compute_angle_residual(th, php, relPos0, relVel, v0, kDrag, P, Fp, mTmp, relTmp, tTmp))
        {
            return false;
        }
        if (!compute_angle_residual(th, phm, relPos0, relVel, v0, kDrag, P, Fm, mTmp, relTmp, tTmp))
        {
            return false;
        }

        Jout[0][1] = (Fp[0] - Fm[0]) / (2.0 * hph);
        Jout[1][1] = wrap_pi(Fp[1] - Fm[1]) / (2.0 * hph);

        return std::isfinite(Jout[0][0]) && std::isfinite(Jout[1][0]) &&
               std::isfinite(Jout[0][1]) && std::isfinite(Jout[1][1]);
    };

    // ----------------------------
    // Initial guess
    // ----------------------------
    double theta, phi;
    initial_guess_vacuum_lead(relPos0, relVel, v0, P.g, theta, phi);
    theta = std::clamp(theta, P.thetaMin, P.thetaMax);
    phi = wrap_pi(phi);

    double F[2];
    double miss = std::numeric_limits<double>::quiet_NaN();
    Vec3 relMissAtStar{};
    double tStar = std::numeric_limits<double>::quiet_NaN();

    if (!compute_angle_residual(theta, phi, relPos0, relVel, v0, kDrag, P, F, miss, relMissAtStar, tStar))
    {
        out.theta = theta;
        out.phi = phi;
        out.miss = miss;
        out.relMissAtStar = relMissAtStar;
        out.tStar = tStar;

        out.report.status = SolveStatus::InitialResidualFailed;
        out.report.message = "InitialResidualFailed: compute_angle_residual returned false.";
        out.report.lastTheta = theta;
        out.report.lastPhi = phi;
        out.report.lastMiss = miss;
        return out;
    }

    // best-so-far
    double thetaBestAll = theta;
    double phiBestAll = phi;
    double missBestAll = miss;
    Vec3 relBestAll = relMissAtStar;
    double tBestAll = tStar;

    // Jacobian init
    double J[2][2];
    if (!jacobian_fd(theta, phi, F, J))
    {
        out.theta = thetaBestAll;
        out.phi = phiBestAll;
        out.miss = missBestAll;
        out.relMissAtStar = relBestAll;
        out.tStar = tBestAll;

        out.report.status = SolveStatus::JacobianFailed;
        out.report.message = "JacobianFailed: initial FD Jacobian evaluation failed.";
        out.report.lastTheta = theta;
        out.report.lastPhi = phi;
        out.report.lastMiss = miss;
        out.report.lastF0 = F[0];
        out.report.lastF1 = F[1];
        return out;
    }

    // ----------------------------
    // Iteration loop
    // ----------------------------
    double lambda = std::clamp(P.lambdaInit, P.lambdaMin, P.lambdaMax);

    for (int it = 0; it < P.maxIter; ++it)
    {
        out.report.iterations = it + 1;

        out.report.lastTheta = theta;
        out.report.lastPhi = phi;
        out.report.lastMiss = miss;
        out.report.lastF0 = F[0];
        out.report.lastF1 = F[1];

        if (miss <= P.tolMiss)
        {
            break;
        }

        const double missOld = miss;
        bool acceptedGlobal = false;

        for (int lt = 0; lt < P.lambdaTries; ++lt)
        {
            out.report.lastLambda = lambda;

            double dtheta, dphi;
            if (!solve_lm_step_2x2(J, F, lambda, dtheta, dphi))
            {
                out.report.status = SolveStatus::LMStepSingular;
                out.report.message = "LMStepSingular: (J'J + lambda I) solve failed (near singular det).";
                goto finalize;
            }

            double alpha = 1.0;

            bool accepted = false;
            bool hasCandidate = false;

            double missBest = miss;
            Vec3 relBest = relMissAtStar;
            double tBest = tStar;
            double thetaBest = theta;
            double phiBest = phi;

            double Fnew[2];

            for (int ls = 0; ls < P.lineSearchTries; ++ls)
            {
                const double thTry = std::clamp(theta + alpha * dtheta, P.thetaMin, P.thetaMax);
                const double phTry = wrap_pi(phi + alpha * dphi);

                out.report.lastAlpha = alpha;

                double mTry;
                Vec3 relTry;
                double tTry;

                if (compute_angle_residual(thTry, phTry, relPos0, relVel, v0, kDrag, P, Fnew, mTry, relTry, tTry))
                {
                    hasCandidate = true;

                    if (!std::isfinite(missBest) || (std::isfinite(mTry) && (mTry < missBest)))
                    {
                        missBest = mTry;
                        relBest = relTry;
                        tBest = tTry;
                        thetaBest = thTry;
                        phiBest = phTry;
                    }

                    if (std::isfinite(mTry) && (mTry <= miss + P.missEps))
                    {
                        accepted = true;
                        break;
                    }
                }
                else
                {
                    // Record cases where residual evaluation fails during line search
                    out.report.status = SolveStatus::ResidualFailedDuringSearch;
                    out.report.message = "ResidualFailedDuringSearch: compute_angle_residual failed during line search.";
                    out.report.lastTheta = thTry;
                    out.report.lastPhi = phTry;
                }

                alpha *= 0.5;
                if (alpha < P.alphaMin)
                {
                    break;
                }
            }

            if (!accepted)
            {
                if (hasCandidate && std::isfinite(missBest) && (missBest < miss))
                {
                    // Re-evaluate and finalize the best candidate
                    if (compute_angle_residual(thetaBest, phiBest, relPos0, relVel, v0, kDrag, P, Fnew, missBest, relBest, tBest))
                    {
                        accepted = true;
                    }
                }
            }

            if (!accepted)
            {
                lambda = std::clamp(lambda * P.lambdaUpMul, P.lambdaMin, P.lambdaMax);
                continue;
            }

            // Step accepted
            out.report.acceptedSteps += 1;

            if (std::isfinite(missBest) && (missBest < missOld - P.missEps))
            {
                lambda = std::clamp(lambda * P.lambdaDownMul, P.lambdaMin, P.lambdaMax);
            }

            // Broyden update
            const double sdx[2] =
            {
                thetaBest - theta,
                wrap_pi(phiBest - phi)
            };

            const double ydf[2] =
            {
                Fnew[0] - F[0],
                wrap_pi(Fnew[1] - F[1])
            };

            broyden_update(J, sdx, ydf, P.broydenMinDenom);

            // State update
            theta = thetaBest;
            phi = wrap_pi(phiBest);

            F[0] = Fnew[0];
            F[1] = Fnew[1];

            miss = missBest;
            relMissAtStar = relBest;
            tStar = tBest;

            // best-so-far update
            if (std::isfinite(miss) && (!std::isfinite(missBestAll) || (miss < missBestAll)))
            {
                thetaBestAll = theta;
                phiBestAll = phi;
                missBestAll = miss;
                relBestAll = relMissAtStar;
                tBestAll = tStar;
            }

            acceptedGlobal = true;
            break;
        }

        if (!acceptedGlobal)
        {
            // All lambdaTries failed
            out.report.status = SolveStatus::LambdaTriesExhausted;
            out.report.message = "LambdaTriesExhausted: no acceptable step within lambdaTries.";
            break;
        }
    }

finalize:
    // Return final (best-so-far) result
    out.theta = thetaBestAll;
    out.phi = phiBestAll;
    out.miss = missBestAll;
    out.relMissAtStar = relBestAll;
    out.tStar = tBestAll;

    // Success check
    out.success = (std::isfinite(out.miss) && (out.miss <= P.tolMiss));

    if (out.success)
    {
        out.report.status = SolveStatus::Ok;
        out.report.message = "Ok";
    }
    else
    {
        // Explicitly set status if it remains Ok
        if (out.report.status == SolveStatus::Ok)
        {
            out.report.status = SolveStatus::MaxIterReached;
            out.report.message = "MaxIterReached: did not satisfy tolMiss within maxIter.";
        }
    }

    // Record final state (for diagnostics)
    out.report.lastTheta = theta;
    out.report.lastPhi = phi;
    out.report.lastMiss = miss;
    out.report.lastLambda = lambda;
    out.report.lastF0 = F[0];
    out.report.lastF1 = F[1];

    return out;
}
