#include <cmath>
#include <stdexcept>
#include <string>
#include <pybind11/pybind11.h>

#include "ballistic_solver_core.hpp"

namespace py = pybind11;

namespace
{
Vec3 to_vec3(const py::sequence& s)
{
    if (py::len(s) != 3)
    {
        throw std::runtime_error("Expected a sequence of length 3.");
    }

    Vec3 v{};
    v.x = py::cast<double>(s[0]);
    v.y = py::cast<double>(s[1]);
    v.z = py::cast<double>(s[2]);
    return v;
}

py::list vec3_to_list(const Vec3& v)
{
    py::list out;
    out.append(v.x);
    out.append(v.y);
    out.append(v.z);
    return out;
}

BallisticParams params_from_dict(const py::dict& d)
{
    BallisticParams p{};
    auto getf = [&](const char* key) { return py::cast<double>(d[key]); };
    auto geti = [&](const char* key) { return py::cast<int>(d[key]); };

    const std::string arc = py::cast<std::string>(d["arcMode"]);
    p.arcMode = (arc == "high") ? ArcMode::High : ArcMode::Low;
    p.g = getf("g");
    p.wind = to_vec3(py::cast<py::sequence>(d["wind"]));
    p.dt = getf("dt");
    p.tMax = getf("tMax");
    p.tolMiss = getf("tolMiss");
    p.beta = getf("beta");
    p.maxIter = geti("maxIter");
    p.lambdaInit = getf("lambdaInit");
    p.lambdaMin = getf("lambdaMin");
    p.lambdaMax = getf("lambdaMax");
    p.lambdaUpMul = getf("lambdaUpMul");
    p.lambdaDownMul = getf("lambdaDownMul");
    p.lambdaTries = geti("lambdaTries");
    p.lineSearchTries = geti("lineSearchTries");
    p.alphaMin = getf("alphaMin");
    p.thetaMin = getf("thetaMin");
    p.thetaMax = getf("thetaMax");
    p.broydenMinDenom = getf("broydenMinDenom");
    p.missEps = getf("missEps");
    p.fdScale = getf("fdScale");
    p.fdMin = getf("fdMin");
    p.fdMax = getf("fdMax");
    return p;
}

void initial_guess_straight_lead_acc(
    const Vec3& relPos0,
    const Vec3& relVel,
    const Vec3& relAcc,
    double v0,
    double& theta0,
    double& phi0)
{
    const double t0 = norm(relPos0) / std::max(v0, 1e-6);
    const Vec3 aim = target_pos_acc(relPos0, relVel, relAcc, t0);
    vec_to_angles(aim, theta0, phi0);
}

bool compute_angle_residual_direct_acc(
    double theta,
    double phi,
    const Vec3& relPos0,
    const Vec3& relVel,
    const Vec3& relAcc,
    double v0,
    double kDrag,
    const BallisticParams& P,
    double F[2],
    double& miss,
    Vec3& relMissAtStarOut,
    double& tStarOut)
{
    const Vec3 projVel0 = v0 * angles_to_dir(theta, phi);
    Vec3 relMissAtStar{};
    double tStar = 0.0;
    find_closest_approach_acc(projVel0, relPos0, relVel, relAcc, kDrag, P, relMissAtStar, tStar);

    miss = norm(relMissAtStar);
    relMissAtStarOut = relMissAtStar;
    tStarOut = tStar;

    const Vec3 aim = target_pos_acc(relPos0, relVel, relAcc, tStar);
    const Vec3 aimCorr = aim - P.beta * relMissAtStar;

    double th0, ph0, th1, ph1;
    vec_to_angles(aim, th0, ph0);
    vec_to_angles(aimCorr, th1, ph1);
    F[0] = th1 - th0;
    F[1] = wrap_pi(ph1 - ph0);
    return std::isfinite(F[0]) && std::isfinite(F[1]) && std::isfinite(miss);
}

double pick_fd_step(const BallisticParams& P, double x)
{
    const double h = P.fdScale * (1.0 + std::fabs(x));
    return std::clamp(h, P.fdMin, P.fdMax);
}

struct CandidateState
{
    double theta = std::numeric_limits<double>::quiet_NaN();
    double phi = std::numeric_limits<double>::quiet_NaN();
    double F[2] = { std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN() };
    double miss = std::numeric_limits<double>::quiet_NaN();
    Vec3 relMissAtStar = { std::numeric_limits<double>::quiet_NaN(),
                           std::numeric_limits<double>::quiet_NaN(),
                           std::numeric_limits<double>::quiet_NaN() };
    double tStar = std::numeric_limits<double>::quiet_NaN();
};

bool eval_candidate_direct(
    CandidateState& out,
    double theta,
    double phi,
    const Vec3& relPos0,
    const Vec3& relVel,
    const Vec3& relAcc,
    double v0,
    double kDrag,
    const BallisticParams& P)
{
    double F[2];
    double miss = std::numeric_limits<double>::quiet_NaN();
    Vec3 relMiss{};
    double tStar = std::numeric_limits<double>::quiet_NaN();
    if (!compute_angle_residual_direct_acc(theta, phi, relPos0, relVel, relAcc, v0, kDrag, P, F, miss, relMiss, tStar))
    {
        return false;
    }

    out.theta = theta;
    out.phi = phi;
    out.F[0] = F[0];
    out.F[1] = F[1];
    out.miss = miss;
    out.relMissAtStar = relMiss;
    out.tStar = tStar;
    return true;
}

bool eval_candidate_aux(
    CandidateState& out,
    double theta,
    double phi,
    const Vec3& relPos0,
    const Vec3& relVel,
    const Vec3& relAcc,
    double v0,
    double kDrag,
    const BallisticParams& P)
{
    double F[2];
    double miss = std::numeric_limits<double>::quiet_NaN();
    Vec3 relMiss{};
    double tStar = std::numeric_limits<double>::quiet_NaN();
    if (!compute_angle_residual_acc(theta, phi, relPos0, relVel, relAcc, v0, kDrag, P, F, miss, relMiss, tStar))
    {
        return false;
    }

    out.theta = theta;
    out.phi = phi;
    out.F[0] = F[0];
    out.F[1] = F[1];
    out.miss = miss;
    out.relMissAtStar = relMiss;
    out.tStar = tStar;
    return true;
}

bool jacobian_direct_fd_acc(
    double theta,
    double phi,
    const double Fbase[2],
    double J[2][2],
    const Vec3& relPos0,
    const Vec3& relVel,
    const Vec3& relAcc,
    double v0,
    double kDrag,
    const BallisticParams& P)
{
    double Fp[2], Fm[2];
    double miss;
    Vec3 relMiss;
    double tStar;

    const double hth = pick_fd_step(P, theta);
    const bool canMinus = theta - hth >= P.thetaMin;
    const bool canPlus = theta + hth <= P.thetaMax;
    if (canMinus && canPlus)
    {
        if (!compute_angle_residual_direct_acc(theta + hth, phi, relPos0, relVel, relAcc, v0, kDrag, P, Fp, miss, relMiss, tStar) ||
            !compute_angle_residual_direct_acc(theta - hth, phi, relPos0, relVel, relAcc, v0, kDrag, P, Fm, miss, relMiss, tStar))
        {
            return false;
        }
        J[0][0] = (Fp[0] - Fm[0]) / (2.0 * hth);
        J[1][0] = wrap_pi(Fp[1] - Fm[1]) / (2.0 * hth);
    }
    else if (canPlus)
    {
        if (!compute_angle_residual_direct_acc(theta + hth, phi, relPos0, relVel, relAcc, v0, kDrag, P, Fp, miss, relMiss, tStar))
        {
            return false;
        }
        J[0][0] = (Fp[0] - Fbase[0]) / hth;
        J[1][0] = wrap_pi(Fp[1] - Fbase[1]) / hth;
    }
    else if (canMinus)
    {
        if (!compute_angle_residual_direct_acc(theta - hth, phi, relPos0, relVel, relAcc, v0, kDrag, P, Fm, miss, relMiss, tStar))
        {
            return false;
        }
        J[0][0] = (Fbase[0] - Fm[0]) / hth;
        J[1][0] = wrap_pi(Fbase[1] - Fm[1]) / hth;
    }
    else
    {
        return false;
    }

    const double hph = pick_fd_step(P, phi);
    if (!compute_angle_residual_direct_acc(theta, wrap_pi(phi + hph), relPos0, relVel, relAcc, v0, kDrag, P, Fp, miss, relMiss, tStar) ||
        !compute_angle_residual_direct_acc(theta, wrap_pi(phi - hph), relPos0, relVel, relAcc, v0, kDrag, P, Fm, miss, relMiss, tStar))
    {
        return false;
    }
    J[0][1] = (Fp[0] - Fm[0]) / (2.0 * hph);
    J[1][1] = wrap_pi(Fp[1] - Fm[1]) / (2.0 * hph);
    return std::isfinite(J[0][0]) && std::isfinite(J[0][1]) && std::isfinite(J[1][0]) && std::isfinite(J[1][1]);
}

py::dict solve_no_aux_impl(
    const Vec3& relPos0,
    const Vec3& relVel,
    const Vec3& relAcc,
    double v0,
    double kDrag,
    const BallisticParams& P)
{
    double theta = 0.0;
    double phi = 0.0;
    initial_guess_straight_lead_acc(relPos0, relVel, relAcc, v0, theta, phi);
    theta = std::clamp(theta, P.thetaMin, P.thetaMax);
    phi = wrap_pi(phi);

    CandidateState cur{};
    if (!eval_candidate_direct(cur, theta, phi, relPos0, relVel, relAcc, v0, kDrag, P))
    {
        throw std::runtime_error("initial direct residual evaluation failed");
    }

    CandidateState best = cur;
    int iterations = 0;
    int acceptedSteps = 0;
    double lambda = std::clamp(P.lambdaInit, P.lambdaMin, P.lambdaMax);

    double J[2][2];
    if (!jacobian_direct_fd_acc(theta, phi, cur.F, J, relPos0, relVel, relAcc, v0, kDrag, P))
    {
        py::dict out;
        out["success"] = std::isfinite(best.miss) && best.miss <= P.tolMiss;
        out["theta"] = best.theta;
        out["phi"] = best.phi;
        out["miss"] = best.miss;
        out["tStar"] = best.tStar;
        out["relMissAtStar"] = vec3_to_list(best.relMissAtStar);
        out["iterations"] = iterations;
        out["acceptedSteps"] = acceptedSteps;
        return out;
    }

    for (int it = 0; it < P.maxIter; ++it)
    {
        iterations = it + 1;
        if (cur.miss <= P.tolMiss)
        {
            break;
        }

        bool accepted = false;
        for (int lt = 0; lt < P.lambdaTries; ++lt)
        {
            double dtheta;
            double dphi;
            if (!solve_lm_step_2x2(J, cur.F, lambda, dtheta, dphi))
            {
                lambda = std::clamp(lambda * P.lambdaUpMul, P.lambdaMin, P.lambdaMax);
                continue;
            }

            double alpha = 1.0;
            CandidateState chosen = cur;
            bool haveChosen = false;
            const double missOld = cur.miss;
            for (int ls = 0; ls < P.lineSearchTries; ++ls)
            {
                const double thetaTry = std::clamp(theta + alpha * dtheta, P.thetaMin, P.thetaMax);
                const double phiTry = wrap_pi(phi + alpha * dphi);
                CandidateState cand{};
                if (eval_candidate_direct(cand, thetaTry, phiTry, relPos0, relVel, relAcc, v0, kDrag, P))
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
                alpha *= 0.5;
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

            acceptedSteps += 1;
            if (chosen.miss < missOld - P.missEps)
            {
                lambda = std::clamp(lambda * P.lambdaDownMul, P.lambdaMin, P.lambdaMax);
            }

            const double sdx[2] = { chosen.theta - theta, wrap_pi(chosen.phi - phi) };
            const double ydf[2] = { chosen.F[0] - cur.F[0], wrap_pi(chosen.F[1] - cur.F[1]) };
            broyden_update(J, sdx, ydf, P.broydenMinDenom);

            theta = chosen.theta;
            phi = wrap_pi(chosen.phi);
            cur = chosen;
            if (cur.miss < best.miss)
            {
                best = cur;
                iterations = it + 1;
            }
            break;
        }

        if (!accepted)
        {
            break;
        }
    }

    py::dict out;
    out["success"] = std::isfinite(best.miss) && best.miss <= P.tolMiss;
    out["theta"] = best.theta;
    out["phi"] = best.phi;
    out["miss"] = best.miss;
    out["tStar"] = best.tStar;
    out["relMissAtStar"] = vec3_to_list(best.relMissAtStar);
    out["iterations"] = iterations;
    out["acceptedSteps"] = acceptedSteps;
    return out;
}

py::dict solve_aux_fixed_point_impl(
    const Vec3& relPos0,
    const Vec3& relVel,
    const Vec3& relAcc,
    double v0,
    double kDrag,
    const BallisticParams& P)
{
    double theta = 0.0;
    double phi = 0.0;
    initial_guess_vacuum_lead_acc(relPos0, relVel, relAcc, v0, P.arcMode, P.g, theta, phi);
    theta = std::clamp(theta, P.thetaMin, P.thetaMax);
    phi = wrap_pi(phi);

    CandidateState cur{};
    if (!eval_candidate_aux(cur, theta, phi, relPos0, relVel, relAcc, v0, kDrag, P))
    {
        throw std::runtime_error("initial auxiliary residual evaluation failed");
    }

    CandidateState best = cur;
    int iterations = 0;
    int acceptedSteps = 0;

    for (int it = 0; it < P.maxIter; ++it)
    {
        iterations = it + 1;
        if (cur.miss <= P.tolMiss)
        {
            break;
        }

        double alpha = 1.0;
        CandidateState chosen = cur;
        bool haveChosen = false;
        bool accepted = false;
        const double missOld = cur.miss;

        for (int ls = 0; ls < P.lineSearchTries; ++ls)
        {
            const double thetaTry = std::clamp(theta + alpha * cur.F[0], P.thetaMin, P.thetaMax);
            const double phiTry = wrap_pi(phi + alpha * cur.F[1]);
            CandidateState cand{};
            if (eval_candidate_aux(cand, thetaTry, phiTry, relPos0, relVel, relAcc, v0, kDrag, P))
            {
                if (!haveChosen || cand.miss < chosen.miss)
                {
                    chosen = cand;
                    haveChosen = true;
                }
                if (cand.miss <= missOld + P.missEps)
                {
                    chosen = cand;
                    accepted = true;
                    break;
                }
            }

            alpha *= 0.5;
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
            break;
        }

        acceptedSteps += 1;
        theta = chosen.theta;
        phi = wrap_pi(chosen.phi);
        cur = chosen;
        if (cur.miss < best.miss)
        {
            best = cur;
            iterations = it + 1;
        }
    }

    py::dict out;
    out["success"] = std::isfinite(best.miss) && best.miss <= P.tolMiss;
    out["theta"] = best.theta;
    out["phi"] = best.phi;
    out["miss"] = best.miss;
    out["tStar"] = best.tStar;
    out["relMissAtStar"] = vec3_to_list(best.relMissAtStar);
    out["iterations"] = iterations;
    out["acceptedSteps"] = acceptedSteps;
    return out;
}
} // namespace

PYBIND11_MODULE(bench_variants_core, m)
{
    m.doc() = "Benchmark-only native solver variants";

    m.def(
        "solve_no_aux",
        [](py::sequence relPos0,
           py::sequence relVel,
           double v0,
           double kDrag,
           py::dict paramsDict,
           py::object relAccObj) -> py::dict
        {
            Vec3 relAcc = { 0.0, 0.0, 0.0 };
            if (!relAccObj.is_none())
            {
                relAcc = to_vec3(py::cast<py::sequence>(relAccObj));
            }
            return solve_no_aux_impl(
                to_vec3(relPos0),
                to_vec3(relVel),
                relAcc,
                v0,
                kDrag,
                params_from_dict(paramsDict));
        },
        py::arg("relPos0"),
        py::arg("relVel"),
        py::arg("v0"),
        py::arg("kDrag"),
        py::arg("params"),
        py::arg("relAcc") = py::none());

    m.def(
        "solve_aux_fixed_point",
        [](py::sequence relPos0,
           py::sequence relVel,
           double v0,
           double kDrag,
           py::dict paramsDict,
           py::object relAccObj) -> py::dict
        {
            Vec3 relAcc = { 0.0, 0.0, 0.0 };
            if (!relAccObj.is_none())
            {
                relAcc = to_vec3(py::cast<py::sequence>(relAccObj));
            }
            return solve_aux_fixed_point_impl(
                to_vec3(relPos0),
                to_vec3(relVel),
                relAcc,
                v0,
                kDrag,
                params_from_dict(paramsDict));
        },
        py::arg("relPos0"),
        py::arg("relVel"),
        py::arg("v0"),
        py::arg("kDrag"),
        py::arg("params"),
        py::arg("relAcc") = py::none());
}
