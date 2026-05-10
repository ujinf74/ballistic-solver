#pragma once

#include <godot_cpp/classes/ref_counted.hpp>
#include <godot_cpp/core/binder_common.hpp>
#include <godot_cpp/core/class_db.hpp>
#include <godot_cpp/variant/dictionary.hpp>
#include <godot_cpp/variant/packed_vector3_array.hpp>
#include <godot_cpp/variant/vector3.hpp>

using namespace godot;

class BallisticSolver : public RefCounted
{
    GDCLASS(BallisticSolver, RefCounted);

protected:
    static void _bind_methods();

public:
    BallisticSolver();
    ~BallisticSolver();

    Dictionary solve(
        const Vector3& rel_pos0,
        const Vector3& rel_vel,
        double v0,
        double k_drag,
        int arc_mode = 0) const;

    PackedVector3Array simulate_from_angles(
        double speed,
        double theta,
        double phi,
        double k_drag,
        double t_max,
        double dt) const;
};
