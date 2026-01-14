#include <cctype>
#include <stdexcept>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

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
                P.arcMode = (am == 0) ? ArcMode::Low : ArcMode::High;
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
}
