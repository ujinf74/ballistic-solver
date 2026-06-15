#pragma once

// Modern C++ surface for the ballistic intercept solver.
//
// This is a source-level C++ API (not the ABI-stable boundary -- that is the C
// API in <ballistic_solver_c_api.h>). Build it into your program with the same
// toolchain as the library. The implementation lives in the compiled library;
// none of the internal headers or helpers leak into your translation unit.

#include <string_view>

#if defined(_WIN32)
    #if defined(BALLISTIC_SOLVER_EXPORTS)
        #define BALLISTIC_SOLVER_API __declspec(dllexport)
    #else
        #define BALLISTIC_SOLVER_API __declspec(dllimport)
    #endif
#else
    #define BALLISTIC_SOLVER_API __attribute__((visibility("default")))
#endif

namespace bs {

struct Vec3
{
    double x{};
    double y{};
    double z{};
};

// What to hit.
struct Problem
{
    Vec3   rel_pos0;     // target position relative to the muzzle at t = 0
    Vec3   rel_vel;      // target velocity relative to the muzzle
    double v0{};         // muzzle speed
    double k_drag{};     // quadratic drag coefficient
    Vec3   rel_acc{};    // optional constant relative acceleration
};

enum class Preset { Fast, Balanced, Precise };
enum class Arc    { Low, High };

// How to solve. Leave the override fields at 0 to use the preset's value.
struct Options
{
    Preset preset   = Preset::Balanced;
    Arc    arc      = Arc::Low;
    Vec3   wind{};                 // air velocity (same frame as the problem)
    double gravity  = 9.80665;

    double tol_miss = 0.0;         // 0 => use the selected preset default
    double dt       = 0.0;         // 0 => use the selected preset default
    double t_max    = 0.0;         // 0 => use the selected preset default
    int    max_iter = 0;           // 0 => use the selected preset default
};

enum class Status
{
    Ok,
    InvalidInput,
    JacobianFailed,
    MaxIterReached,
    NumericalFailure,
};

// The result. A best-effort solution is always returned; check `success`.
struct Intercept
{
    bool             success{};
    double           theta{};                 // elevation, radians
    double           phi{};                   // azimuth, radians
    double           miss{};                  // closest-approach distance
    double           t_star{};                // time of closest approach
    Vec3             rel_miss_at_star{};
    Status           status = Status::InvalidInput;
    std::string_view message{};               // static-lifetime diagnostic (no allocation)
    int              iterations{};

    explicit operator bool() const noexcept { return success; }
};

// Default coordinate-residual solver.
BALLISTIC_SOLVER_API Intercept solve(const Problem& problem, const Options& options = {});

// Auxiliary-residual (ASIR) variant; kept for compatibility/reproducibility.
BALLISTIC_SOLVER_API Intercept solve_aux(const Problem& problem, const Options& options = {});

} // namespace bs
