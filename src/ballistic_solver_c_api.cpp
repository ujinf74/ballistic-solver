#include <cmath>
#include <cstdio>
#include <cstring>

#include "ballistic_solver_c_api.h"
#include "ballistic_solver_core.hpp"

static inline bool finite3_ptr(const double* a)
{
    return a &&
           std::isfinite(a[0]) &&
           std::isfinite(a[1]) &&
           std::isfinite(a[2]);
}

static inline bool finite3_arr(const double a[3])
{
    return finite3_ptr(a);
}

void ballistic_inputs_init(BallisticInputs* in)
{
    if (in == nullptr)
    {
        return;
    }

    in->arcMode = 0;
    in->g       = 9.80665;
    in->dt      = 0.01;
    in->tMax    = 20.0;
    in->tolMiss = 1e-2;
    in->maxIter = 20;
}

int32_t ballistic_solve(const BallisticInputs* in, BallisticOutputs* out)
{
    if (in == nullptr || out == nullptr)
    {
        if (out != nullptr)
        {
            std::snprintf(out->message, sizeof(out->message), "Null pointer argument");
        }
        return -1;
    }

    std::memset(out, 0, sizeof(*out));
    std::snprintf(out->message, sizeof(out->message), "Unknown");

    try
    {
        BallisticParams P;

        if (in->arcMode == 1)
        {
            P.arcMode = ArcMode::High;
        }
        else
        {
            P.arcMode = ArcMode::Low;
        }

        if (std::isfinite(in->g) && in->g > 0.0)       { P.g = in->g; }
        if (std::isfinite(in->dt) && in->dt > 0.0)     { P.dt = in->dt; }
        if (std::isfinite(in->tMax) && in->tMax > 0.0) { P.tMax = in->tMax; }
        if (std::isfinite(in->tolMiss) && in->tolMiss > 0.0) { P.tolMiss = in->tolMiss; }
        if (in->maxIter > 0) { P.maxIter = static_cast<int>(in->maxIter); }

        const Vec3 relPos0 = { in->relPos0[0], in->relPos0[1], in->relPos0[2] };
        const Vec3 relVel  = { in->relVel[0],  in->relVel[1],  in->relVel[2]  };
        P.wind = { in->wind[0], in->wind[1], in->wind[2] };

        const SolverResult r = solve_launch_angles(relPos0, relVel, in->v0, in->kDrag, P);

        out->success = r.success ? 1 : 0;
        out->status  = static_cast<int32_t>(r.report.status);

        out->theta = r.theta;
        out->phi   = r.phi;
        out->miss  = r.miss;
        out->tStar = r.tStar;

        out->relMissAtStar[0] = r.relMissAtStar.x;
        out->relMissAtStar[1] = r.relMissAtStar.y;
        out->relMissAtStar[2] = r.relMissAtStar.z;

        std::snprintf(out->message, sizeof(out->message), "%s", r.report.message.c_str());

        return 0;
    }
    catch (...)
    {
        out->success = 0;
        out->status = static_cast<int32_t>(SolveStatus::InvalidInput);
        std::snprintf(out->message, sizeof(out->message), "Exception caught in ballistic_solve");
        return -2;
    }
}

uint32_t ballistic_solver_abi_version(void)
{
    return (uint32_t)BALLISTIC_SOLVER_ABI_VERSION;
}

const char* ballistic_solver_version_string(void)
{
    return BALLISTIC_SOLVER_VERSION_STRING;
}

void ballistic_rk4_step(
    double r3[3],
    double v3[3],
    double h,
    double kDrag,
    double g,
    const double wind3[3])
{
    if (!r3 || !v3 || !wind3) { return; }
    if (!finite3_arr(r3) || !finite3_arr(v3) || !finite3_arr(wind3)) { return; }
    if (!std::isfinite(h) || h <= 0.0) { return; }
    if (!std::isfinite(kDrag)) { return; }
    if (!std::isfinite(g) || g <= 0.0) { return; }

    try
    {
        State s{};
        s.r = { r3[0], r3[1], r3[2] };
        s.v = { v3[0], v3[1], v3[2] };

        const Vec3 wind = { wind3[0], wind3[1], wind3[2] };

        rk4_step(s, h, kDrag, g, wind);

        r3[0] = s.r.x; r3[1] = s.r.y; r3[2] = s.r.z;
        v3[0] = s.v.x; v3[1] = s.v.y; v3[2] = s.v.z;
    }
    catch (...)
    {
        // C ABI: swallow exceptions
    }
}

int32_t ballistic_simulate_trajectory(
    const double r0_3[3],
    const double v0_3[3],
    double kDrag,
    double g,
    const double wind3[3],
    double dt,
    int32_t steps,
    double* out_pos3,
    double* out_vel3,
    double* out_t,
    int32_t* out_count)
{
    if (!r0_3 || !v0_3 || !wind3 || !out_pos3) { return -1; }
    if (!finite3_arr(r0_3) || !finite3_arr(v0_3) || !finite3_arr(wind3)) { return -2; }
    if (!std::isfinite(kDrag)) { return -3; }
    if (!std::isfinite(g) || g <= 0.0) { return -4; }
    if (!std::isfinite(dt) || dt <= 0.0) { return -5; }
    if (steps < 0) { return -6; }

    try
    {
        State s{};
        s.r = { r0_3[0], r0_3[1], r0_3[2] };
        s.v = { v0_3[0], v0_3[1], v0_3[2] };

        const Vec3 wind = { wind3[0], wind3[1], wind3[2] };

        const int32_t count = steps + 1;
        if (out_count) { *out_count = count; }

        double t = 0.0;

        // sample 0
        out_pos3[0] = s.r.x; out_pos3[1] = s.r.y; out_pos3[2] = s.r.z;
        if (out_vel3)
        {
            out_vel3[0] = s.v.x; out_vel3[1] = s.v.y; out_vel3[2] = s.v.z;
        }
        if (out_t) { out_t[0] = t; }

        for (int32_t i = 0; i < steps; ++i)
        {
            rk4_step(s, dt, kDrag, g, wind);
            t += dt;

            const int32_t j = (i + 1) * 3;
            out_pos3[j + 0] = s.r.x;
            out_pos3[j + 1] = s.r.y;
            out_pos3[j + 2] = s.r.z;

            if (out_vel3)
            {
                out_vel3[j + 0] = s.v.x;
                out_vel3[j + 1] = s.v.y;
                out_vel3[j + 2] = s.v.z;
            }
            if (out_t) { out_t[i + 1] = t; }
        }

        return 0;
    }
    catch (...)
    {
        return -7;
    }
}

int32_t ballistic_simulate_trajectory_from_angles(
    const double r0_3[3],
    double speed,
    double theta,
    double phi,
    double kDrag,
    double g,
    const double wind3[3],
    double dt,
    int32_t steps,
    double* out_pos3,
    double* out_vel3,
    double* out_t,
    int32_t* out_count)
{
    if (!r0_3 || !wind3 || !out_pos3) { return -1; }
    if (!finite3_arr(r0_3) || !finite3_arr(wind3)) { return -2; }
    if (!std::isfinite(speed) || speed <= 0.0) { return -3; }
    if (!std::isfinite(theta) || !std::isfinite(phi)) { return -4; }

    const Vec3 dir = angles_to_dir(theta, phi);
    const double v0_3[3] = { speed * dir.x, speed * dir.y, speed * dir.z };

    return ballistic_simulate_trajectory(
        r0_3, v0_3,
        kDrag, g, wind3,
        dt, steps,
        out_pos3, out_vel3, out_t, out_count);
}

int32_t ballistic_find_closest_approach(
    const BallisticInputs* in,
    double theta,
    double phi,
    double* out_tStar,
    double out_relMissAtStar3[3],
    double* out_miss)
{
    if (!in) { return -1; }
    if (!out_tStar || !out_relMissAtStar3 || !out_miss) { return -2; }

    if (!std::isfinite(theta) || !std::isfinite(phi)) { return -3; }
    if (!std::isfinite(in->v0) || in->v0 <= 0.0) { return -4; }
    if (!std::isfinite(in->kDrag)) { return -5; }
    if (!std::isfinite(in->g) || in->g <= 0.0) { return -6; }
    if (!std::isfinite(in->dt) || in->dt <= 0.0) { return -7; }
    if (!std::isfinite(in->tMax) || in->tMax <= 0.0) { return -8; }
    if (!finite3_ptr(in->relPos0) || !finite3_ptr(in->relVel) || !finite3_ptr(in->wind)) { return -9; }

    try
    {
        BallisticParams P;
        P.arcMode = (in->arcMode == 1) ? ArcMode::High : ArcMode::Low;
        P.g = in->g;
        P.dt = in->dt;
        P.tMax = in->tMax;
        P.wind = { in->wind[0], in->wind[1], in->wind[2] };

        const Vec3 relPos0 = { in->relPos0[0], in->relPos0[1], in->relPos0[2] };
        const Vec3 relVel  = { in->relVel[0],  in->relVel[1],  in->relVel[2]  };

        const Vec3 dir = angles_to_dir(theta, phi);
        const Vec3 projVel0 = in->v0 * dir;

        Vec3 relMissAtStar{};
        double tStar = 0.0;

        find_closest_approach(projVel0, relPos0, relVel, in->kDrag, P, relMissAtStar, tStar);

        out_relMissAtStar3[0] = relMissAtStar.x;
        out_relMissAtStar3[1] = relMissAtStar.y;
        out_relMissAtStar3[2] = relMissAtStar.z;

        *out_tStar = tStar;
        *out_miss = norm(relMissAtStar);

        return 0;
    }
    catch (...)
    {
        return -10;
    }
}

int32_t ballistic_vacuum_arc_angles_to_point(
    const double R3[3],
    double v0,
    int32_t arcMode,
    double g,
    double* out_theta,
    double* out_phi)
{
    if (!R3 || !out_theta || !out_phi) { return -1; }
    if (!finite3_arr(R3)) { return -2; }
    if (!std::isfinite(v0) || v0 <= 0.0) { return -3; }
    if (!std::isfinite(g) || g <= 0.0) { return -4; }

    try
    {
        const Vec3 R = { R3[0], R3[1], R3[2] };
        const ArcMode mode = (arcMode == 1) ? ArcMode::High : ArcMode::Low;

        double th = 0.0, ph = 0.0;
        const bool ok = vacuum_arc_angles_to_point(R, v0, mode, g, th, ph);

        *out_theta = th;
        *out_phi = ph;

        return ok ? 1 : 0;
    }
    catch (...)
    {
        return -5;
    }
}

void ballistic_initial_guess_vacuum_lead(
    const double relPos0_3[3],
    const double relVel_3[3],
    double v0,
    int32_t arcMode,
    double g,
    double* out_theta,
    double* out_phi)
{
    if (!relPos0_3 || !relVel_3 || !out_theta || !out_phi) { return; }
    if (!finite3_arr(relPos0_3) || !finite3_arr(relVel_3)) { return; }
    if (!std::isfinite(v0) || v0 <= 0.0) { return; }
    if (!std::isfinite(g) || g <= 0.0) { return; }

    try
    {
        const Vec3 relPos0 = { relPos0_3[0], relPos0_3[1], relPos0_3[2] };
        const Vec3 relVel  = { relVel_3[0],  relVel_3[1],  relVel_3[2]  };

        const ArcMode mode = (arcMode == 1) ? ArcMode::High : ArcMode::Low;

        double th = 0.0, ph = 0.0;
        initial_guess_vacuum_lead(relPos0, relVel, v0, mode, g, th, ph);

        *out_theta = th;
        *out_phi = ph;
    }
    catch (...)
    {
        // C ABI: swallow exceptions
    }
}
