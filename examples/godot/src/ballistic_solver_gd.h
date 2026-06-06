#pragma once

#include <godot_cpp/classes/ref_counted.hpp>
#include <godot_cpp/core/binder_common.hpp>
#include <godot_cpp/core/class_db.hpp>
#include <godot_cpp/variant/dictionary.hpp>
#include <godot_cpp/variant/packed_vector3_array.hpp>
#include <godot_cpp/variant/vector3.hpp>

#include "ballistic_solver_core.hpp"

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

// Position-only constant-acceleration tracker, wrapping the C++ core directly
// (no C ABI). Ingests (timestamp, relative position) measurements and produces
// lead predictions / firing solutions via the predictor seam.
class BallisticTracker : public RefCounted
{
    GDCLASS(BallisticTracker, RefCounted);

    TargetTracker tracker_;

protected:
    static void _bind_methods();

public:
    BallisticTracker();
    ~BallisticTracker();

    // Reset the filter and set process (jerk) and measurement (position) noise.
    void configure(double process_noise, double meas_noise);

    void update(double t, const Vector3& rel_pos);

    Vector3 predict(double tau) const;
    Vector3 estimated_position() const;
    Vector3 estimated_velocity() const;
    Vector3 estimated_acceleration() const;

    // Lead-fire firing solution against the current track (predictor seam).
    Dictionary solve(double v0, double k_drag, int arc_mode = 0) const;
};
