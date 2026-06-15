// Implementation of the modern C++ surface declared in <ballistic_solver.hpp>.
// Maps the public bs:: types onto the internal (global-namespace) header-only
// core and back, so consumers never see the core's symbols.

#include "ballistic_solver.hpp"

#include "ballistic_solver_core.hpp"   // BallisticParams, SolverResult, solve_launch_angles, ...
#include "bs/solve_coord_lead.hpp"     // solve_launch_angles_coord_lead (default core)

namespace {

::BallisticParams to_params(const bs::Options& o)
{
    const ::ParamPreset pp =
        o.preset == bs::Preset::Fast    ? ::ParamPreset::Fast :
        o.preset == bs::Preset::Precise ? ::ParamPreset::Precise :
                                          ::ParamPreset::Balanced;

    ::BallisticParams P = ::make_params_preset(pp);
    P.arcMode = (o.arc == bs::Arc::High) ? ::ArcMode::High : ::ArcMode::Low;
    P.g = o.gravity;
    P.wind = ::Vec3{ o.wind.x, o.wind.y, o.wind.z };

    if (o.tol_miss > 0.0) P.tolMiss = o.tol_miss;   // 0 => keep preset value
    if (o.dt > 0.0)       P.dt = o.dt;
    if (o.t_max > 0.0)    P.tMax = o.t_max;
    if (o.max_iter > 0)   P.maxIter = o.max_iter;
    return P;
}

bs::Status map_status(::SolveStatus s)
{
    switch (s)
    {
    case ::SolveStatus::Ok:             return bs::Status::Ok;
    case ::SolveStatus::InvalidInput:   return bs::Status::InvalidInput;
    case ::SolveStatus::JacobianFailed: return bs::Status::JacobianFailed;
    case ::SolveStatus::MaxIterReached: return bs::Status::MaxIterReached;
    default:                            return bs::Status::NumericalFailure;  // solve_aux LM paths
    }
}

std::string_view message_for(bs::Status s)
{
    switch (s)
    {
    case bs::Status::Ok:               return "ok";
    case bs::Status::InvalidInput:     return "invalid input (check v0, gravity, dt, t_max, max_iter, arc range)";
    case bs::Status::JacobianFailed:   return "vacuum Jacobian could not be formed at the seed";
    case bs::Status::MaxIterReached:   return "did not reach tol_miss within max_iter (best-effort result returned)";
    case bs::Status::NumericalFailure: return "numerical failure during the auxiliary-residual solve";
    }
    return "";
}

bs::Intercept run(const bs::Problem& p, const bs::Options& o, bool aux)
{
    const ::BallisticParams P = to_params(o);
    const ::Vec3 rp{ p.rel_pos0.x, p.rel_pos0.y, p.rel_pos0.z };
    const ::Vec3 rv{ p.rel_vel.x, p.rel_vel.y, p.rel_vel.z };
    const ::Vec3 ra{ p.rel_acc.x, p.rel_acc.y, p.rel_acc.z };

    const ::SolverResult res = aux
        ? ::solve_launch_angles(rp, rv, p.v0, p.k_drag, P, ra)
        : ::solve_launch_angles_coord_lead(rp, rv, p.v0, p.k_drag, P, ra);

    bs::Intercept out;
    out.success          = res.success;
    out.theta            = res.theta;
    out.phi              = res.phi;
    out.miss             = res.miss;
    out.t_star           = res.tStar;
    out.rel_miss_at_star = { res.relMissAtStar.x, res.relMissAtStar.y, res.relMissAtStar.z };
    out.status           = map_status(res.report.status);
    out.message          = message_for(out.status);
    out.iterations       = res.report.iterations;
    return out;
}

} // namespace

namespace bs {

Intercept solve(const Problem& problem, const Options& options)
{
    return run(problem, options, /*aux=*/false);
}

Intercept solve_aux(const Problem& problem, const Options& options)
{
    return run(problem, options, /*aux=*/true);
}

} // namespace bs
