#include <cctype>
#include <stdexcept>
#include <string>
#include <pybind11/pybind11.h>

#include "ballistic_solver_core.hpp"

namespace py = pybind11;

static Vec3 to_vec3(const py::sequence& s)
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

static py::list vec3_to_list(const Vec3& v)
{
    py::list out;
    out.append(v.x);
    out.append(v.y);
    out.append(v.z);
    return out;
}

static int normalize_arc_mode(py::object arcMode)
{
    if (py::isinstance<py::str>(arcMode))
    {
        std::string s = py::cast<std::string>(arcMode);
        for (char& c : s)
        {
            c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
        }

        if (s == "low" || s == "l" || s == "0")
        {
            return 0;
        }
        if (s == "high" || s == "h" || s == "1")
        {
            return 1;
        }

        throw std::runtime_error("Invalid arcMode string. Use 'low' or 'high'.");
    }

    return py::cast<int>(arcMode);
}

PYBIND11_MODULE(_core, m)
{
    m.doc() = "Ballistic solver (pybind11 core)";

    py::enum_<ArcMode>(m, "ArcMode")
        .value("Low", ArcMode::Low)
        .value("High", ArcMode::High);

    py::class_<BallisticParams>(m, "BallisticParams")
        .def(py::init<>())

        .def_readwrite("arcMode", &BallisticParams::arcMode)
        .def_readwrite("g", &BallisticParams::g)
        .def_property(
            "wind",
            [](const BallisticParams& p) -> py::list
            {
                return vec3_to_list(p.wind);
            },
            [](BallisticParams& p, const py::sequence& wind)
            {
                p.wind = to_vec3(wind);
            })
        .def_readwrite("dt", &BallisticParams::dt)
        .def_readwrite("tMax", &BallisticParams::tMax)
        .def_readwrite("tolMiss", &BallisticParams::tolMiss)
        .def_readwrite("beta", &BallisticParams::beta)
        .def_readwrite("maxIter", &BallisticParams::maxIter)
        .def_readwrite("lambdaInit", &BallisticParams::lambdaInit)
        .def_readwrite("lambdaMin", &BallisticParams::lambdaMin)
        .def_readwrite("lambdaMax", &BallisticParams::lambdaMax)
        .def_readwrite("lambdaUpMul", &BallisticParams::lambdaUpMul)
        .def_readwrite("lambdaDownMul", &BallisticParams::lambdaDownMul)
        .def_readwrite("lambdaTries", &BallisticParams::lambdaTries)
        .def_readwrite("lineSearchTries", &BallisticParams::lineSearchTries)
        .def_readwrite("alphaMin", &BallisticParams::alphaMin)
        .def_readwrite("thetaMin", &BallisticParams::thetaMin)
        .def_readwrite("thetaMax", &BallisticParams::thetaMax)
        .def_readwrite("broydenMinDenom", &BallisticParams::broydenMinDenom)
        .def_readwrite("missEps", &BallisticParams::missEps)
        .def_readwrite("fdScale", &BallisticParams::fdScale)
        .def_readwrite("fdMin", &BallisticParams::fdMin)
        .def_readwrite("fdMax", &BallisticParams::fdMax)
        .def_readwrite("gsMaxIter", &BallisticParams::gsMaxIter)
        .def_readwrite("gsTolAbs", &BallisticParams::gsTolAbs)
        .def_readwrite("gsTolRel", &BallisticParams::gsTolRel);

    m.def(
        "solve",
        [](py::sequence relPos0,
           py::sequence relVel,
           double v0,
           double kDrag,
           py::object arcMode,
           py::object paramsObj) -> py::dict
        {
            BallisticParams P{};

            if (!paramsObj.is_none())
            {
                P = py::cast<BallisticParams>(paramsObj);
            }

            if (!arcMode.is_none())
            {
                const int am = normalize_arc_mode(arcMode);
                P.arcMode = (am == 1) ? ArcMode::High : ArcMode::Low;
            }

            const Vec3 r0 = to_vec3(relPos0);
            const Vec3 rv = to_vec3(relVel);

            SolverResult r = solve_launch_angles(r0, rv, v0, kDrag, P);

            py::dict out;
            out["success"] = r.success;
            out["theta"] = r.theta;
            out["phi"] = r.phi;
            out["miss"] = r.miss;
            out["tStar"] = r.tStar;
            out["relMissAtStar"] = vec3_to_list(r.relMissAtStar);

            out["status"] = static_cast<int>(r.report.status);
            out["message"] = r.report.message;
            out["iterations"] = r.report.iterations;
            out["acceptedSteps"] = r.report.acceptedSteps;
            out["lastLambda"] = r.report.lastLambda;
            out["lastAlpha"] = r.report.lastAlpha;

            return out;
        },
        py::arg("relPos0"),
        py::arg("relVel"),
        py::arg("v0"),
        py::arg("kDrag"),
        py::arg("arcMode") = py::none(),
        py::arg("params") = py::none());

    m.def(
        "rk4_step",
        [](py::sequence r,
           py::sequence v,
           double h,
           double kDrag,
           py::object paramsObj) -> py::dict
        {
            BallisticParams P{};
            if (!paramsObj.is_none()) P = py::cast<BallisticParams>(paramsObj);

            State s{};
            s.r = to_vec3(r);
            s.v = to_vec3(v);

            rk4_step(s, h, kDrag, P.g, P.wind);

            py::dict out;
            out["r"] = vec3_to_list(s.r);
            out["v"] = vec3_to_list(s.v);
            return out;
        },
        py::arg("r"),
        py::arg("v"),
        py::arg("h"),
        py::arg("kDrag"),
        py::arg("params") = py::none());

    m.def(
        "simulate",
        [](py::sequence r0,
           py::sequence v0,
           double kDrag,
           int steps,
           py::object paramsObj) -> py::dict
        {
            if (steps < 0) throw std::runtime_error("steps must be >= 0.");

            BallisticParams P{};
            if (!paramsObj.is_none()) P = py::cast<BallisticParams>(paramsObj);

            State s{};
            s.r = to_vec3(r0);
            s.v = to_vec3(v0);

            py::list ts, rs, vs;
            double t = 0.0;

            ts.append(t);
            rs.append(vec3_to_list(s.r));
            vs.append(vec3_to_list(s.v));

            for (int i = 0; i < steps; ++i)
            {
                rk4_step(s, P.dt, kDrag, P.g, P.wind);
                t += P.dt;
                ts.append(t);
                rs.append(vec3_to_list(s.r));
                vs.append(vec3_to_list(s.v));
            }

            py::dict out;
            out["t"] = ts;
            out["r"] = rs;
            out["v"] = vs;
            return out;
        },
        py::arg("r0"),
        py::arg("v0"),
        py::arg("kDrag"),
        py::arg("steps"),
        py::arg("params") = py::none());

    m.def(
        "closest_approach",
        [](py::sequence relPos0,
           py::sequence relVel,
           double v0,
           double kDrag,
           double theta,
           double phi,
           py::object paramsObj) -> py::dict
        {
            BallisticParams P{};
            if (!paramsObj.is_none()) P = py::cast<BallisticParams>(paramsObj);

            const Vec3 dir = angles_to_dir(theta, phi);
            const Vec3 projVel0 = v0 * dir;

            Vec3 relMiss{};
            double tStar = 0.0;

            find_closest_approach(projVel0, to_vec3(relPos0), to_vec3(relVel), kDrag, P, relMiss, tStar);

            py::dict out;
            out["tStar"] = tStar;
            out["miss"] = norm(relMiss);
            out["relMissAtStar"] = vec3_to_list(relMiss);
            return out;
        },
        py::arg("relPos0"),
        py::arg("relVel"),
        py::arg("v0"),
        py::arg("kDrag"),
        py::arg("theta"),
        py::arg("phi"),
        py::arg("params") = py::none());

    m.def(
        "vacuum_arc_angles_to_point",
        [](py::sequence R,
           double v0,
           py::object arcMode,
           double g) -> py::dict
        {
            const int am = normalize_arc_mode(arcMode);
            const ArcMode mode = (am == 1) ? ArcMode::High : ArcMode::Low;

            double theta = 0.0, phi = 0.0;
            const bool ok = vacuum_arc_angles_to_point(to_vec3(R), v0, mode, g, theta, phi);

            py::dict out;
            out["reachable"] = ok;
            out["theta"] = theta;
            out["phi"] = phi;
            return out;
        },
        py::arg("R"),
        py::arg("v0"),
        py::arg("arcMode") = py::str("low"),
        py::arg("g") = 9.80665);

    m.def(
        "vacuum_lead_initial_guess",
        [](py::sequence relPos0,
           py::sequence relVel,
           double v0,
           py::object arcMode,
           double g) -> py::dict
        {
            const int am = normalize_arc_mode(arcMode);
            const ArcMode mode = (am == 1) ? ArcMode::High : ArcMode::Low;

            double theta = 0.0, phi = 0.0;
            initial_guess_vacuum_lead(to_vec3(relPos0), to_vec3(relVel), v0, mode, g, theta, phi);

            py::dict out;
            out["theta"] = theta;
            out["phi"] = phi;
            return out;
        },
        py::arg("relPos0"),
        py::arg("relVel"),
        py::arg("v0"),
        py::arg("arcMode") = py::str("low"),
        py::arg("g") = 9.80665);
}
