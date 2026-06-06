#include "ballistic_solver_gd.h"

#include <gdextension_interface.h>

#include <cmath>
#include <vector>

#include <godot_cpp/core/class_db.hpp>
#include <godot_cpp/core/defs.hpp>
#include <godot_cpp/godot.hpp>

#include "ballistic_solver_c_api.h"

using namespace godot;

namespace
{
void godot_to_solver(const Vector3& src, double dst[3])
{
    dst[0] = src.x;
    dst[1] = src.z;
    dst[2] = src.y;
}

Vector3 solver_to_godot(const double src[3])
{
    return Vector3(
        static_cast<real_t>(src[0]),
        static_cast<real_t>(src[2]),
        static_cast<real_t>(src[1]));
}

Vec3 godot_to_vec3(const Vector3& v)
{
    double a[3];
    godot_to_solver(v, a);
    return Vec3{ a[0], a[1], a[2] };
}

Vector3 vec3_to_godot(const Vec3& v)
{
    const double a[3] = { v.x, v.y, v.z };
    return solver_to_godot(a);
}
}

BallisticSolver::BallisticSolver() = default;

BallisticSolver::~BallisticSolver() = default;

void BallisticSolver::_bind_methods()
{
    ClassDB::bind_method(
        D_METHOD("solve", "rel_pos0", "rel_vel", "v0", "k_drag", "arc_mode"),
        &BallisticSolver::solve,
        DEFVAL(0));
    ClassDB::bind_method(
        D_METHOD("simulate_from_angles", "speed", "theta", "phi", "k_drag", "t_max", "dt"),
        &BallisticSolver::simulate_from_angles);
}

Dictionary BallisticSolver::solve(
    const Vector3& rel_pos0,
    const Vector3& rel_vel,
    double v0,
    double k_drag,
    int arc_mode) const
{
    BallisticInputs in;
    BallisticOutputs out{};
    ballistic_inputs_init(&in);

    godot_to_solver(rel_pos0, in.relPos0);
    godot_to_solver(rel_vel, in.relVel);
    in.v0 = v0;
    in.kDrag = k_drag;
    in.arcMode = arc_mode;

    const int rc = ballistic_solve(&in, &out);

    Dictionary d;
    d["call_rc"] = rc;
    if (rc < 0)
    {
        d["success"] = false;
        d["status"] = -1;
        d["message"] = "ballistic_solve API call failed";
        return d;
    }

    d["success"] = out.success != 0;
    d["status"] = out.status;
    d["theta"] = out.theta;
    d["phi"] = out.phi;
    d["miss"] = out.miss;
    d["t_star"] = out.tStar;
    d["rel_miss_at_star"] = solver_to_godot(out.relMissAtStar);
    d["iterations"] = out.iterations;
    d["accepted_steps"] = out.acceptedSteps;
    d["last_lambda"] = out.lastLambda;
    d["last_alpha"] = out.lastAlpha;
    d["message"] = String::utf8(out.message);
    return d;
}

PackedVector3Array BallisticSolver::simulate_from_angles(
    double speed,
    double theta,
    double phi,
    double k_drag,
    double t_max,
    double dt) const
{
    PackedVector3Array points;
    if (speed <= 0.0 || t_max <= 0.0 || dt <= 0.0)
    {
        return points;
    }

    const int steps = static_cast<int>(std::ceil(t_max / dt));
    std::vector<double> pos(static_cast<size_t>(steps + 1) * 3);
    std::vector<double> vel(static_cast<size_t>(steps + 1) * 3);
    std::vector<double> times(static_cast<size_t>(steps + 1));
    const double r0[3] = {0.0, 0.0, 0.0};
    const double wind[3] = {0.0, 0.0, 0.0};
    int32_t out_count = 0;

    const int rc = ballistic_simulate_trajectory_from_angles(
        r0,
        speed,
        theta,
        phi,
        k_drag,
        9.80665,
        wind,
        dt,
        steps,
        pos.data(),
        vel.data(),
        times.data(),
        &out_count);

    if (rc != 0)
    {
        return points;
    }

    points.resize(out_count);
    for (int i = 0; i < out_count; ++i)
    {
        const int j = 3 * i;
        points.set(i, solver_to_godot(&pos[j]));
    }
    return points;
}

BallisticTracker::BallisticTracker() = default;

BallisticTracker::~BallisticTracker() = default;

void BallisticTracker::_bind_methods()
{
    ClassDB::bind_method(D_METHOD("configure", "process_noise", "meas_noise"), &BallisticTracker::configure);
    ClassDB::bind_method(D_METHOD("update", "t", "rel_pos"), &BallisticTracker::update);
    ClassDB::bind_method(D_METHOD("predict", "tau"), &BallisticTracker::predict);
    ClassDB::bind_method(D_METHOD("estimated_position"), &BallisticTracker::estimated_position);
    ClassDB::bind_method(D_METHOD("estimated_velocity"), &BallisticTracker::estimated_velocity);
    ClassDB::bind_method(D_METHOD("estimated_acceleration"), &BallisticTracker::estimated_acceleration);
    ClassDB::bind_method(
        D_METHOD("solve", "v0", "k_drag", "arc_mode"),
        &BallisticTracker::solve,
        DEFVAL(0));
}

void BallisticTracker::configure(double process_noise, double meas_noise)
{
    tracker_ = TargetTracker(process_noise, meas_noise);
}

void BallisticTracker::update(double t, const Vector3& rel_pos)
{
    tracker_.update(t, godot_to_vec3(rel_pos));
}

Vector3 BallisticTracker::predict(double tau) const
{
    return vec3_to_godot(tracker_.predict(tau));
}

Vector3 BallisticTracker::estimated_position() const
{
    return vec3_to_godot(tracker_.position());
}

Vector3 BallisticTracker::estimated_velocity() const
{
    return vec3_to_godot(tracker_.velocity());
}

Vector3 BallisticTracker::estimated_acceleration() const
{
    return vec3_to_godot(tracker_.acceleration());
}

Dictionary BallisticTracker::solve(double v0, double k_drag, int arc_mode) const
{
    BallisticParams P;
    P.arcMode = (arc_mode == 1) ? ArcMode::High : ArcMode::Low;

    const SolverResult r = solve_launch_angles_predicted(tracker_.predictor(), v0, k_drag, P);

    Dictionary d;
    d["success"] = r.success;
    d["status"] = static_cast<int>(r.report.status);
    d["theta"] = r.theta;
    d["phi"] = r.phi;
    d["miss"] = r.miss;
    d["t_star"] = r.tStar;
    d["rel_miss_at_star"] = Vector3(
        static_cast<real_t>(r.relMissAtStar.x),
        static_cast<real_t>(r.relMissAtStar.z),
        static_cast<real_t>(r.relMissAtStar.y));
    d["iterations"] = r.report.iterations;
    return d;
}

void initialize_ballistic_solver_module(ModuleInitializationLevel level)
{
    if (level != MODULE_INITIALIZATION_LEVEL_SCENE)
    {
        return;
    }

    ClassDB::register_class<BallisticSolver>();
    ClassDB::register_class<BallisticTracker>();
}

void uninitialize_ballistic_solver_module(ModuleInitializationLevel level)
{
    if (level != MODULE_INITIALIZATION_LEVEL_SCENE)
    {
        return;
    }
}

extern "C"
{
GDExtensionBool GDE_EXPORT ballistic_solver_gdextension_init(
    GDExtensionInterfaceGetProcAddress get_proc_address,
    GDExtensionClassLibraryPtr library,
    GDExtensionInitialization* initialization)
{
    GDExtensionBinding::InitObject init_obj(get_proc_address, library, initialization);
    init_obj.register_initializer(initialize_ballistic_solver_module);
    init_obj.register_terminator(uninitialize_ballistic_solver_module);
    init_obj.set_minimum_library_initialization_level(MODULE_INITIALIZATION_LEVEL_SCENE);
    return init_obj.init();
}
}
