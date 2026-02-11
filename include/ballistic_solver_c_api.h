#pragma once

#include <stdint.h>

#if defined(_WIN32)
    #if defined(BALLISTIC_SOLVER_EXPORTS)
        #define BALLISTIC_SOLVER_C_API __declspec(dllexport)
    #else
        #define BALLISTIC_SOLVER_C_API __declspec(dllimport)
    #endif
    #define BALLISTIC_SOLVER_CALL __cdecl
#else
    #define BALLISTIC_SOLVER_C_API __attribute__((visibility("default")))
    #define BALLISTIC_SOLVER_CALL
#endif

#ifdef __cplusplus
extern "C" {
#endif

#ifndef BALLISTIC_SOLVER_ABI_VERSION
    #define BALLISTIC_SOLVER_ABI_VERSION 2u
#endif

#ifndef BALLISTIC_SOLVER_VERSION_STRING
    #define BALLISTIC_SOLVER_VERSION_STRING "unknown"
#endif

BALLISTIC_SOLVER_C_API uint32_t BALLISTIC_SOLVER_CALL ballistic_solver_abi_version(void);
BALLISTIC_SOLVER_C_API const char* BALLISTIC_SOLVER_CALL ballistic_solver_version_string(void);

typedef struct BallisticInputs
{
    double relPos0[3];
    double relVel[3];
    double v0;
    double kDrag;

    int32_t arcMode; // 0=Low, 1=High
    int32_t _pad0;
    double g;
    double wind[3];
    double dt;
    double tMax;
    double tolMiss;
    int32_t maxIter;
    int32_t _pad1;
} BallisticInputs;

typedef struct BallisticOutputs
{
    int32_t success;
    int32_t status;

    double theta;
    double phi;
    double miss;
    double tStar;

    double relMissAtStar[3];

    char message[256];
} BallisticOutputs;

BALLISTIC_SOLVER_C_API void ballistic_inputs_init(BallisticInputs* in);

BALLISTIC_SOLVER_C_API int32_t ballistic_solve(const BallisticInputs* in, BallisticOutputs* out);

BALLISTIC_SOLVER_C_API void BALLISTIC_SOLVER_CALL ballistic_rk4_step(
    double r3[3],
    double v3[3],
    double h,
    double kDrag,
    double g,
    const double wind3[3]);

BALLISTIC_SOLVER_C_API int32_t BALLISTIC_SOLVER_CALL ballistic_simulate_trajectory(
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
    int32_t* out_count);

BALLISTIC_SOLVER_C_API int32_t BALLISTIC_SOLVER_CALL ballistic_simulate_trajectory_from_angles(
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
    int32_t* out_count);

BALLISTIC_SOLVER_C_API int32_t BALLISTIC_SOLVER_CALL ballistic_find_closest_approach(
    const BallisticInputs* in,
    double theta,
    double phi,
    double* out_tStar,
    double out_relMissAtStar3[3],
    double* out_miss);

BALLISTIC_SOLVER_C_API int32_t BALLISTIC_SOLVER_CALL ballistic_vacuum_arc_angles_to_point(
    const double R3[3],
    double v0,
    int32_t arcMode, /* 0=Low, 1=High */
    double g,
    double* out_theta,
    double* out_phi);

BALLISTIC_SOLVER_C_API void BALLISTIC_SOLVER_CALL ballistic_initial_guess_vacuum_lead(
    const double relPos0_3[3],
    const double relVel_3[3],
    double v0,
    int32_t arcMode, /* 0=Low, 1=High */
    double g,
    double* out_theta,
    double* out_phi);

#ifdef __cplusplus
}
#endif