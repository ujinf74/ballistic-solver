#include "ballistic_solver_gd.h"

#include <gdextension_interface.h>

#include <godot_cpp/core/class_db.hpp>
#include <godot_cpp/core/defs.hpp>
#include <godot_cpp/godot.hpp>

#include "ballistic_solver_c_api.h"

using namespace godot;

BallisticSolver::BallisticSolver() = default;

BallisticSolver::~BallisticSolver() = default;

void BallisticSolver::_bind_methods()
{
    ClassDB::bind_method(
        D_METHOD("solve", "rel_pos0", "rel_vel", "v0", "k_drag", "arc_mode"),
        &BallisticSolver::solve,
        DEFVAL(0));
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

    in.relPos0[0] = rel_pos0.x;
    in.relPos0[1] = rel_pos0.y;
    in.relPos0[2] = rel_pos0.z;
    in.relVel[0] = rel_vel.x;
    in.relVel[1] = rel_vel.y;
    in.relVel[2] = rel_vel.z;
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
    d["rel_miss_at_star"] = Vector3(
        static_cast<real_t>(out.relMissAtStar[0]),
        static_cast<real_t>(out.relMissAtStar[1]),
        static_cast<real_t>(out.relMissAtStar[2]));
    d["iterations"] = out.iterations;
    d["accepted_steps"] = out.acceptedSteps;
    d["last_lambda"] = out.lastLambda;
    d["last_alpha"] = out.lastAlpha;
    d["message"] = String::utf8(out.message);
    return d;
}

void initialize_ballistic_solver_module(ModuleInitializationLevel level)
{
    if (level != MODULE_INITIALIZATION_LEVEL_SCENE)
    {
        return;
    }

    ClassDB::register_class<BallisticSolver>();
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
