#pragma once

#include "residual.hpp"

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
// Step helpers (LM + line search + lambda tries)
// ================================================================
struct CandidateState
{
    double theta = std::numeric_limits<double>::quiet_NaN();
    double phi = std::numeric_limits<double>::quiet_NaN();

    double F[2] = { std::numeric_limits<double>::quiet_NaN(),
                    std::numeric_limits<double>::quiet_NaN() };

    double miss = std::numeric_limits<double>::quiet_NaN();
    Vec3 relMissAtStar = { std::numeric_limits<double>::quiet_NaN(),
                           std::numeric_limits<double>::quiet_NaN(),
                           std::numeric_limits<double>::quiet_NaN() };
    double tStar = std::numeric_limits<double>::quiet_NaN();
};

struct StepResult
{
    bool accepted = false;
    bool hadCandidate = false;

    CandidateState best;

    double dtheta = std::numeric_limits<double>::quiet_NaN();
    double dphi = std::numeric_limits<double>::quiet_NaN();

    double lastAlpha = std::numeric_limits<double>::quiet_NaN();
};

inline void init_candidate_from_current(
    CandidateState& c,
    double theta,
    double phi,
    const double F[2],
    double miss,
    const Vec3& relMissAtStar,
    double tStar)
{
    c.theta = theta;
    c.phi = phi;
    c.F[0] = F[0];
    c.F[1] = F[1];
    c.miss = miss;
    c.relMissAtStar = relMissAtStar;
    c.tStar = tStar;
}

inline bool evaluate_candidate(
    CandidateState& out,
    double thetaTry,
    double phiTry,
    const Vec3& relPos0,
    const Vec3& relVel,
    const Vec3& relAcc,
    double v0,
    double kDrag,
    const BallisticParams& P)
{
    double Fnew[2];
    double mTry;
    Vec3 relTry;
    double tTry;

    if (!compute_angle_residual_acc(thetaTry, phiTry, relPos0, relVel, relAcc, v0, kDrag, P, Fnew, mTry, relTry, tTry))
    {
        return false;
    }

    out.theta = thetaTry;
    out.phi = phiTry;

    out.F[0] = Fnew[0];
    out.F[1] = Fnew[1];

    out.miss = mTry;
    out.relMissAtStar = relTry;
    out.tStar = tTry;

    return std::isfinite(out.miss) && std::isfinite(out.F[0]) && std::isfinite(out.F[1]);
}

inline StepResult line_search_best(
    double theta,
    double phi,
    const double F[2],
    double miss,
    const Vec3& relMissAtStar,
    double tStar,
    double dtheta,
    double dphi,
    const Vec3& relPos0,
    const Vec3& relVel,
    const Vec3& relAcc,
    double v0,
    double kDrag,
    const BallisticParams& P,
    SolveReport& report)
{
    StepResult res{};
    res.dtheta = dtheta;
    res.dphi = dphi;

    init_candidate_from_current(res.best, theta, phi, F, miss, relMissAtStar, tStar);

    double alpha = 1.0;

    for (int ls = 0; ls < P.lineSearchTries; ++ls)
    {
        const double thTry = std::clamp(theta + alpha * dtheta, P.thetaMin, P.thetaMax);
        const double phTry = wrap_pi(phi + alpha * dphi);

        report.lastAlpha = alpha;
        res.lastAlpha = alpha;

        CandidateState cand{};
        if (evaluate_candidate(cand, thTry, phTry, relPos0, relVel, relAcc, v0, kDrag, P))
        {
            res.hadCandidate = true;

            if (!std::isfinite(res.best.miss) || (std::isfinite(cand.miss) && (cand.miss < res.best.miss)))
            {
                res.best = cand;
            }

            if (std::isfinite(cand.miss) && (cand.miss <= miss + P.missEps))
            {
                res.accepted = true;
                return res;
            }
        }
        else
        {
            report.status = SolveStatus::ResidualFailedDuringSearch;
            report.message = "ResidualFailedDuringSearch: compute_angle_residual failed during line search.";
            report.lastTheta = thTry;
            report.lastPhi = phTry;
        }

        alpha *= P.lineSearchShrink;
        if (alpha < P.alphaMin)
        {
            break;
        }
    }

    if (!res.accepted && res.hadCandidate && std::isfinite(res.best.miss) && (res.best.miss < miss))
    {
        CandidateState reeval{};
        if (evaluate_candidate(reeval, res.best.theta, res.best.phi, relPos0, relVel, relAcc, v0, kDrag, P))
        {
            res.best = reeval;
            res.accepted = true;
        }
    }

    return res;
}

inline bool accept_step(
    double& theta,
    double& phi,
    double F[2],
    double& miss,
    Vec3& relMissAtStar,
    double& tStar,
    double J[2][2],
    double& lambda,
    const BallisticParams& P,
    const StepResult& step,
    double missOld,
    SolveReport& report)
{
    if (!step.accepted)
    {
        return false;
    }

    report.acceptedSteps += 1;

    if (std::isfinite(step.best.miss) && (step.best.miss < missOld - P.missEps))
    {
        lambda = std::clamp(lambda * P.lambdaDownMul, P.lambdaMin, P.lambdaMax);
    }

    const double sdx[2] =
    {
        step.best.theta - theta,
        wrap_pi(step.best.phi - phi)
    };

    const double ydf[2] =
    {
        step.best.F[0] - F[0],
        wrap_pi(step.best.F[1] - F[1])
    };

    broyden_update(J, sdx, ydf, P.broydenMinDenom);

    theta = step.best.theta;
    phi = wrap_pi(step.best.phi);

    F[0] = step.best.F[0];
    F[1] = step.best.F[1];

    miss = step.best.miss;
    relMissAtStar = step.best.relMissAtStar;
    tStar = step.best.tStar;

    return true;
}

inline bool try_lm_step(
    double& theta,
    double& phi,
    double F[2],
    double& miss,
    Vec3& relMissAtStar,
    double& tStar,
    double J[2][2],
    double& lambda,
    const Vec3& relPos0,
    const Vec3& relVel,
    const Vec3& relAcc,
    double v0,
    double kDrag,
    const BallisticParams& P,
    SolveReport& report)
{
    report.lastLambda = lambda;

    double dtheta;
    double dphi;

    if (!solve_lm_step_2x2(J, F, lambda, dtheta, dphi))
    {
        report.status = SolveStatus::LMStepSingular;
        report.message = "LMStepSingular: (J'J + lambda I) solve failed (near singular det).";
        return false;
    }

    const double missOld = miss;

    StepResult step = line_search_best(
        theta, phi, F, miss, relMissAtStar, tStar,
        dtheta, dphi,
        relPos0, relVel, relAcc, v0, kDrag, P,
        report);

    if (!step.accepted)
    {
        lambda = std::clamp(lambda * P.lambdaUpMul, P.lambdaMin, P.lambdaMax);
        return false;
    }

    return accept_step(
        theta, phi, F, miss, relMissAtStar, tStar,
        J, lambda, P,
        step, missOld,
        report);
}

inline Vec3 vacuum_required_displacement(
    const Vec3& relPos0,
    const Vec3& relVel,
    const Vec3& relAcc,
    double g,
    double t)
{
    return target_pos_acc(relPos0, relVel, relAcc, t) + Vec3{ 0.0, 0.0, 0.5 * g * t * t };
}

template <std::size_t N>
inline void add_time_grid_seeds(
    std::array<CandidateState, N>& seeds,
    int& seedCount,
    const Vec3& relPos0,
    const Vec3& relVel,
    const Vec3& relAcc,
    double v0,
    double kDrag,
    const BallisticParams& P)
{
    for (int i = 3; i <= 5 && seedCount < static_cast<int>(N); ++i)
    {
        const double u = static_cast<double>(i) / 10.0;
        const double t = std::max(P.dt, P.tMax * (0.2 + 0.8 * u));
        const Vec3 aim = vacuum_required_displacement(relPos0, relVel, relAcc, P.g, t);

        double theta;
        double phi;
        vec_to_angles(aim, theta, phi);
        theta = std::clamp(theta, P.thetaMin, P.thetaMax);
        phi = wrap_pi(phi);

        CandidateState cand{};
        if (evaluate_candidate(cand, theta, phi, relPos0, relVel, relAcc, v0, kDrag, P))
        {
            seeds[seedCount++] = cand;
        }
    }
}

inline bool solve_auxiliary_multistart(
    double& bestTheta,
    double& bestPhi,
    double& bestMiss,
    Vec3& bestRelMiss,
    double& bestTime,
    const Vec3& relPos0,
    const Vec3& relVel,
    const Vec3& relAcc,
    double v0,
    double kDrag,
    const BallisticParams& P,
    SolveReport& report)
{
    double theta = 0.0;
    double phi = 0.0;
    initial_guess_vacuum_lead_acc(relPos0, relVel, relAcc, v0, P.arcMode, P.g, theta, phi, P.tMax);
    theta = std::clamp(theta, P.thetaMin, P.thetaMax);
    phi = wrap_pi(phi);

    std::array<CandidateState, 4> seeds{};
    int seedCount = 0;

    CandidateState cur{};
    if (evaluate_candidate(cur, theta, phi, relPos0, relVel, relAcc, v0, kDrag, P))
    {
        seeds[seedCount++] = cur;
    }

    add_time_grid_seeds(seeds, seedCount, relPos0, relVel, relAcc, v0, kDrag, P);

    if (seedCount == 0)
    {
        return false;
    }

    std::sort(seeds.begin(), seeds.begin() + seedCount, [](const CandidateState& a, const CandidateState& b)
    {
        return a.miss < b.miss;
    });

    auto pick_fd_step = [&](double x) -> double
    {
        const double h = P.fdScale * (1.0 + std::fabs(x));
        return std::clamp(h, P.fdMin, P.fdMax);
    };

    auto jacobian_aux_fd = [&](double th, double ph, const double Fbase[2], double Jout[2][2]) -> bool
    {
        double Fp[2], Fm[2];
        double mTmp;
        Vec3 relTmp;
        double tTmp;

        const double hth = pick_fd_step(th);
        const bool canPlus = (th + hth <= P.thetaMax);
        const bool canMinus = (th - hth >= P.thetaMin);

        if (canMinus && canPlus)
        {
            if (!compute_angle_residual_acc(th + hth, ph, relPos0, relVel, relAcc, v0, kDrag, P, Fp, mTmp, relTmp, tTmp) ||
                !compute_angle_residual_acc(th - hth, ph, relPos0, relVel, relAcc, v0, kDrag, P, Fm, mTmp, relTmp, tTmp))
            {
                return false;
            }
            Jout[0][0] = (Fp[0] - Fm[0]) / (2.0 * hth);
            Jout[1][0] = wrap_pi(Fp[1] - Fm[1]) / (2.0 * hth);
        }
        else if (canPlus)
        {
            if (!compute_angle_residual_acc(th + hth, ph, relPos0, relVel, relAcc, v0, kDrag, P, Fp, mTmp, relTmp, tTmp))
            {
                return false;
            }
            Jout[0][0] = (Fp[0] - Fbase[0]) / hth;
            Jout[1][0] = wrap_pi(Fp[1] - Fbase[1]) / hth;
        }
        else if (canMinus)
        {
            if (!compute_angle_residual_acc(th - hth, ph, relPos0, relVel, relAcc, v0, kDrag, P, Fm, mTmp, relTmp, tTmp))
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
        if (!compute_angle_residual_acc(th, wrap_pi(ph + hph), relPos0, relVel, relAcc, v0, kDrag, P, Fp, mTmp, relTmp, tTmp) ||
            !compute_angle_residual_acc(th, wrap_pi(ph - hph), relPos0, relVel, relAcc, v0, kDrag, P, Fm, mTmp, relTmp, tTmp))
        {
            return false;
        }

        Jout[0][1] = (Fp[0] - Fm[0]) / (2.0 * hph);
        Jout[1][1] = wrap_pi(Fp[1] - Fm[1]) / (2.0 * hph);
        return std::isfinite(Jout[0][0]) && std::isfinite(Jout[1][0]) &&
               std::isfinite(Jout[0][1]) && std::isfinite(Jout[1][1]);
    };

    CandidateState best = seeds[0];
    bool ranAnySeed = false;

    auto run_seed = [&](const CandidateState& seed) -> void
    {
        CandidateState current = seed;
        CandidateState localBest = current;
        double thetaSeed = current.theta;
        double phiSeed = current.phi;
        double lambda = std::clamp(P.lambdaInit, P.lambdaMin, P.lambdaMax);

        double J[2][2];
        if (!jacobian_aux_fd(thetaSeed, phiSeed, current.F, J))
        {
            return;
        }

        ranAnySeed = true;

        for (int it = 0; it < P.maxIter; ++it)
        {
            if (current.miss <= P.tolMiss)
            {
                break;
            }

            bool accepted = false;
            for (int lt = 0; lt < P.lambdaTries; ++lt)
            {
                double dtheta;
                double dphi;
                if (!solve_lm_step_2x2(J, current.F, lambda, dtheta, dphi))
                {
                    lambda = std::clamp(lambda * P.lambdaUpMul, P.lambdaMin, P.lambdaMax);
                    continue;
                }

                double alpha = 1.0;
                CandidateState chosen = current;
                bool haveChosen = false;
                const double missOld = current.miss;
                for (int ls = 0; ls < P.lineSearchTries; ++ls)
                {
                    const double thetaTry = std::clamp(thetaSeed + alpha * dtheta, P.thetaMin, P.thetaMax);
                    const double phiTry = wrap_pi(phiSeed + alpha * dphi);
                    CandidateState cand{};
                    if (evaluate_candidate(cand, thetaTry, phiTry, relPos0, relVel, relAcc, v0, kDrag, P))
                    {
                        if (!haveChosen || cand.miss < chosen.miss)
                        {
                            chosen = cand;
                            haveChosen = true;
                        }
                        if (cand.miss <= missOld + P.missEps)
                        {
                            accepted = true;
                            chosen = cand;
                            break;
                        }
                    }
                    alpha *= P.lineSearchShrink;
                    if (alpha < P.alphaMin)
                    {
                        break;
                    }
                }

                if (!accepted && haveChosen && chosen.miss < missOld)
                {
                    accepted = true;
                }
                if (!accepted)
                {
                    lambda = std::clamp(lambda * P.lambdaUpMul, P.lambdaMin, P.lambdaMax);
                    continue;
                }

                report.acceptedSteps += 1;
                if (chosen.miss < missOld - P.missEps)
                {
                    lambda = std::clamp(lambda * P.lambdaDownMul, P.lambdaMin, P.lambdaMax);
                }

                const double sdx[2] = { chosen.theta - thetaSeed, wrap_pi(chosen.phi - phiSeed) };
                const double ydf[2] = { chosen.F[0] - current.F[0], wrap_pi(chosen.F[1] - current.F[1]) };
                broyden_update(J, sdx, ydf, P.broydenMinDenom);

                thetaSeed = chosen.theta;
                phiSeed = wrap_pi(chosen.phi);
                current = chosen;
                if (current.miss < localBest.miss)
                {
                    localBest = current;
                }
                break;
            }

            if (!accepted)
            {
                break;
            }
        }

        if (localBest.miss < best.miss)
        {
            best = localBest;
        }
    };

    for (int i = 0; i < seedCount; ++i)
    {
        run_seed(seeds[i]);
        if (best.miss <= P.tolMiss)
        {
            break;
        }
    }

    if (ranAnySeed && std::isfinite(best.miss) && best.miss < bestMiss)
    {
        bestTheta = best.theta;
        bestPhi = best.phi;
        bestMiss = best.miss;
        bestRelMiss = best.relMissAtStar;
        bestTime = best.tStar;
        return true;
    }
    return false;
}
