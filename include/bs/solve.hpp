#pragma once

#include "lm.hpp"

// ================================================================
// Solve launch angles (Broyden)
// ================================================================
inline SolverResult solve_launch_angles(
    const Vec3& relPos0,
    const Vec3& relVel,
    double v0,
    double kDrag,
    const BallisticParams& P = BallisticParams{},
    const Vec3& relAcc = Vec3{ 0.0, 0.0, 0.0 })
{
    SolverResult out{};

    // ----------------------------
    // Input validation
    // ----------------------------
    if (!std::isfinite(v0) || v0 <= 0.0 || !std::isfinite(kDrag) ||
        !std::isfinite(P.g) || P.g <= 0.0 || !std::isfinite(P.dt) || P.dt <= 0.0 ||
        !std::isfinite(P.tMax) || P.tMax <= 0.0 || P.maxIter <= 0 ||
        !std::isfinite(P.lineSearchShrink) || P.lineSearchShrink <= 0.0 || P.lineSearchShrink >= 1.0 ||
        !std::isfinite(P.beta) || P.beta <= 0.0 ||
        P.thetaMin >= P.thetaMax)
    {
        out.report.status = SolveStatus::InvalidInput;
        out.report.message = "InvalidInput: v0/g/dt/tMax/maxIter/theta range check failed.";
        return out;
    }

    BallisticParams residualP = P;

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
            if (!compute_angle_residual_acc(th + hth, ph, relPos0, relVel, relAcc, v0, kDrag, residualP, Fp, mTmp, relTmp, tTmp))
            {
                return false;
            }
            if (!compute_angle_residual_acc(th - hth, ph, relPos0, relVel, relAcc, v0, kDrag, residualP, Fm, mTmp, relTmp, tTmp))
            {
                return false;
            }

            Jout[0][0] = (Fp[0] - Fm[0]) / (2.0 * hth);
            Jout[1][0] = wrap_pi(Fp[1] - Fm[1]) / (2.0 * hth);
        }
        else if (canPlus)
        {
            if (!compute_angle_residual_acc(th + hth, ph, relPos0, relVel, relAcc, v0, kDrag, residualP, Fp, mTmp, relTmp, tTmp))
            {
                return false;
            }

            Jout[0][0] = (Fp[0] - Fbase[0]) / hth;
            Jout[1][0] = wrap_pi(Fp[1] - Fbase[1]) / hth;
        }
        else if (canMinus)
        {
            if (!compute_angle_residual_acc(th - hth, ph, relPos0, relVel, relAcc, v0, kDrag, residualP, Fm, mTmp, relTmp, tTmp))
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

        if (!compute_angle_residual_acc(th, php, relPos0, relVel, relAcc, v0, kDrag, residualP, Fp, mTmp, relTmp, tTmp))
        {
            return false;
        }
        if (!compute_angle_residual_acc(th, phm, relPos0, relVel, relAcc, v0, kDrag, residualP, Fm, mTmp, relTmp, tTmp))
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
    initial_guess_vacuum_lead_acc(relPos0, relVel, relAcc, v0, P.arcMode, P.g, theta, phi, P.tMax);
    theta = std::clamp(theta, P.thetaMin, P.thetaMax);
    phi = wrap_pi(phi);

    double F[2];
    double miss = std::numeric_limits<double>::quiet_NaN();
    Vec3 relMissAtStar{};
    double tStar = std::numeric_limits<double>::quiet_NaN();

    if (!compute_angle_residual_acc(theta, phi, relPos0, relVel, relAcc, v0, kDrag, residualP, F, miss, relMissAtStar, tStar))
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
    double bestTheta = theta;
    double bestPhi = phi;
    double bestMiss = miss;
    Vec3 bestRelMiss = relMissAtStar;
    double bestTime = tStar;

    if (miss > P.tolMiss)
    {
        CandidateState preStep{};
        double preDtheta;
        double preDphi;
        const Vec3 preAim = target_pos_acc(relPos0, relVel, relAcc, tStar);
        const bool havePreStep = compute_auxiliary_delta(
            preAim, relMissAtStar, P.preStepBeta, v0, P, preDtheta, preDphi);
        const double thetaTry = std::clamp(theta + (havePreStep ? preDtheta : F[0]), P.thetaMin, P.thetaMax);
        const double phiTry = wrap_pi(phi + (havePreStep ? preDphi : F[1]));
        if (evaluate_candidate(preStep, thetaTry, phiTry, relPos0, relVel, relAcc, v0, kDrag, residualP) &&
            std::isfinite(preStep.miss) && preStep.miss < miss)
        {
            theta = preStep.theta;
            phi = wrap_pi(preStep.phi);

            F[0] = preStep.F[0];
            F[1] = preStep.F[1];

            miss = preStep.miss;
            relMissAtStar = preStep.relMissAtStar;
            tStar = preStep.tStar;

            bestTheta = theta;
            bestPhi = phi;
            bestMiss = miss;
            bestRelMiss = relMissAtStar;
            bestTime = tStar;
        }
    }

    double J[2][2];
    if (!jacobian_fd(theta, phi, F, J))
    {
        out.theta = bestTheta;
        out.phi = bestPhi;
        out.miss = bestMiss;
        out.relMissAtStar = bestRelMiss;
        out.tStar = bestTime;

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

    bool hardFailure = false;

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

        bool acceptedGlobal = false;
        for (int lt = 0; lt < P.lambdaTries; ++lt)
        {
            if (try_lm_step(
                    theta, phi, F, miss, relMissAtStar, tStar,
                    J, lambda,
                    relPos0, relVel, relAcc, v0, kDrag, residualP,
                    out.report))
            {
                acceptedGlobal = true;

                if (std::isfinite(miss) && (!std::isfinite(bestMiss) || (miss < bestMiss)))
                {
                    bestTheta = theta;
                    bestPhi = phi;
                    bestMiss = miss;
                    bestRelMiss = relMissAtStar;
                    bestTime = tStar;
                }

                break;
            }

            if (out.report.status == SolveStatus::LMStepSingular)
            {
                hardFailure = true;
                break;
            }
        }

        if (hardFailure)
        {
            break;
        }

        if (!acceptedGlobal)
        {
            double auxDtheta;
            double auxDphi;
            const Vec3 auxAim = target_pos_acc(relPos0, relVel, relAcc, tStar);
            if (compute_auxiliary_delta(auxAim, relMissAtStar, P.preStepBeta, v0, P, auxDtheta, auxDphi))
            {
                const double missOld = miss;
                StepResult auxStep = line_search_best(
                    theta, phi, F, miss, relMissAtStar, tStar,
                    auxDtheta, auxDphi,
                    relPos0, relVel, relAcc, v0, kDrag, residualP,
                    out.report);
                acceptedGlobal = accept_step(
                    theta, phi, F, miss, relMissAtStar, tStar,
                    J, lambda, residualP,
                    auxStep, missOld,
                    out.report);
                if (acceptedGlobal && std::isfinite(miss) && (!std::isfinite(bestMiss) || (miss < bestMiss)))
                {
                    bestTheta = theta;
                    bestPhi = phi;
                    bestMiss = miss;
                    bestRelMiss = relMissAtStar;
                    bestTime = tStar;
                }
            }

            if (!acceptedGlobal)
            {
                out.report.status = SolveStatus::LambdaTriesExhausted;
                out.report.message = "LambdaTriesExhausted: no acceptable step within lambdaTries.";
                break;
            }
        }
    }

    if (!std::isfinite(bestMiss) || bestMiss > P.tolMiss)
    {
        solve_auxiliary_multistart(
            bestTheta, bestPhi, bestMiss, bestRelMiss, bestTime,
            relPos0, relVel, relAcc, v0, kDrag, residualP,
            out.report);
    }

    // Return final (best-so-far) result
    out.theta = bestTheta;
    out.phi = bestPhi;
    out.miss = bestMiss;
    out.relMissAtStar = bestRelMiss;
    out.tStar = bestTime;

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
