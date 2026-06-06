#pragma once

#include "math_utils.hpp"

// ================================================================
// // Parameters (configurable)
// ================================================================
struct BallisticParams
{
    ArcMode arcMode = ArcMode::Low; // Arc mode
    
    double g = 9.80665;        // Gravity acceleration
    Vec3 wind = { 0.0, 0.0, 0.0 };        // Wind vector

    double dt = 0.01;          // RK4 basic step
    double tMax = 40.0;        // Maximum simulation time

    double tolMiss = 1e-2;     // Success tolerance (miss distance)
    double beta = 0.04;        // Residual target offset correction factor
    double preStepBeta = 1.0;  // One-shot pre-LM correction factor

    int maxIter = 20;

    double lambdaInit = 1e-6;
    double lambdaMin = 1e-12;
    double lambdaMax = 1e+6;
    double lambdaUpMul = 10.0;
    double lambdaDownMul = 0.3;
    int lambdaTries = 2;

    int lineSearchTries = 1;
    double lineSearchShrink = 0.5;
    double alphaMin = 1e-4;

    double thetaMin = -M_PI / 2.0;
    double thetaMax =  M_PI / 2.0;

    double broydenMinDenom = 1e-12;

    double missEps = 1e-12;

    double fdScale = 1e-4;
    double fdMin = 1e-6;
    double fdMax = 1e-3;

    int gsMaxIter = 20;
    double gsTolAbs = 1e-8;
    double gsTolRel = 1e-8;
};

inline BallisticParams make_params_preset(ParamPreset preset)
{
    BallisticParams P{};

    switch (preset)
    {
    case ParamPreset::Fast:
        P.dt = 0.02;
        P.tolMiss = 5e-2;
        P.maxIter = 12;
        P.lambdaTries = 4;
        P.lineSearchTries = 8;
        P.gsMaxIter = 12;
        P.gsTolAbs = 1e-7;
        P.gsTolRel = 1e-7;
        break;

    case ParamPreset::Precise:
        P.dt = 0.01;
        P.tolMiss = 1e-5;
        P.maxIter = 80;
        P.lambdaTries = 10;
        P.lineSearchTries = 16;
        P.gsMaxIter = 40;
        P.gsTolAbs = 1e-10;
        P.gsTolRel = 1e-10;
        P.fdScale = 1e-5;
        P.fdMin = 1e-7;
        P.fdMax = 1e-4;
        break;

    case ParamPreset::Balanced:
    default:
        break;
    }

    return P;
}

inline bool k_drag_from_physical(
    double airDensity,
    double dragCoefficient,
    double area,
    double mass,
    double& kDrag)
{
    if (!std::isfinite(airDensity) || !std::isfinite(dragCoefficient) ||
        !std::isfinite(area) || !std::isfinite(mass) ||
        airDensity < 0.0 || dragCoefficient < 0.0 || area < 0.0 || mass <= 0.0)
    {
        kDrag = std::numeric_limits<double>::quiet_NaN();
        return false;
    }

    kDrag = 0.5 * airDensity * dragCoefficient * area / mass;
    return std::isfinite(kDrag);
}

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
