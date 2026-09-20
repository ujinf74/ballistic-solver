#pragma once

#include "lm.hpp"

// Guards every field solve_launch_angles dereferences before it integrates.
inline bool solve_inputs_are_valid(double v0, double kDrag, const BallisticParams& P)
{
    return std::isfinite(v0) && v0 > 0.0 && std::isfinite(kDrag) &&
           std::isfinite(P.g) && P.g > 0.0 && std::isfinite(P.dt) && P.dt > 0.0 &&
           std::isfinite(P.tMax) && P.tMax > 0.0 && P.maxIter > 0 &&
           std::isfinite(P.lineSearchShrink) && P.lineSearchShrink > 0.0 && P.lineSearchShrink < 1.0 &&
           std::isfinite(P.beta) && P.beta > 0.0 &&
           P.thetaMin < P.thetaMax;
}

// One auxiliary-residual correction taken before the LM loop starts. Applied
// (and reported) only when it strictly improves the miss.
inline bool try_auxiliary_prestep(
    double& theta,
    double& phi,
    double F[2],
    double& miss,
    Vec3& relMissAtStar,
    double& tStar,
    const Vec3& relPos0,
    const Vec3& relVel,
    const Vec3& relAcc,
    double v0,
    double kDrag,
    const BallisticParams& P)
{
    CandidateState preStep{};
    double preDtheta;
    double preDphi;
    const Vec3 preAim = target_pos_acc(relPos0, relVel, relAcc, tStar);
    const bool havePreStep = compute_auxiliary_delta(
        preAim, relMissAtStar, P.preStepBeta, v0, P, preDtheta, preDphi);
    const double thetaTry = std::clamp(theta + (havePreStep ? preDtheta : F[0]), P.thetaMin, P.thetaMax);
    const double phiTry = wrap_pi(phi + (havePreStep ? preDphi : F[1]));

    if (evaluate_candidate(preStep, thetaTry, phiTry, relPos0, relVel, relAcc, v0, kDrag, P) &&
        std::isfinite(preStep.miss) && preStep.miss < miss)
    {
        theta = preStep.theta;
        phi = wrap_pi(preStep.phi);

        F[0] = preStep.F[0];
        F[1] = preStep.F[1];

        miss = preStep.miss;
        relMissAtStar = preStep.relMissAtStar;
        tStar = preStep.tStar;
        return true;
    }

    return false;
}

// ================================================================
// Solve launch angles (Broyden)
// ================================================================
inline SolverResult solve_launch_angles(
    const Vec3& relPos0,
    const Vec3& relVel,
    double v0,
    double kDrag,
    const BallisticParams& P = BallisticParams{},
    const Vec3& relAcc = Vec3{ 0.0, 0.0, 0.0 },
    bool allowMultistart = true)
{
    SolverResult out{};

    if (!solve_inputs_are_valid(v0, kDrag, P))
    {
        out.report.status = SolveStatus::InvalidInput;
        out.report.message = "InvalidInput: v0/g/dt/tMax/maxIter/theta range check failed.";
        return out;
    }

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

    if (!compute_angle_residual_acc(theta, phi, relPos0, relVel, relAcc, v0, kDrag, P, F, miss, relMissAtStar, tStar))
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

    if (miss > P.tolMiss &&
        try_auxiliary_prestep(theta, phi, F, miss, relMissAtStar, tStar,
            relPos0, relVel, relAcc, v0, kDrag, P))
    {
        bestTheta = theta;
        bestPhi = phi;
        bestMiss = miss;
        bestRelMiss = relMissAtStar;
        bestTime = tStar;
    }

    double J[2][2];
    if (!jacobian_angles_fd(theta, phi, F, relPos0, relVel, relAcc, v0, kDrag, P, J))
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
                    relPos0, relVel, relAcc, v0, kDrag, P,
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
                    relPos0, relVel, relAcc, v0, kDrag, P,
                    out.report);
                acceptedGlobal = accept_step(
                    theta, phi, F, miss, relMissAtStar, tStar,
                    J, lambda, P,
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

    if (allowMultistart && (!std::isfinite(bestMiss) || bestMiss > P.tolMiss))
    {
        solve_auxiliary_multistart(
            bestTheta, bestPhi, bestMiss, bestRelMiss, bestTime,
            relPos0, relVel, relAcc, v0, kDrag, P,
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
