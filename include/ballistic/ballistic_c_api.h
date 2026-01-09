#pragma once

#include <stdint.h>

#if defined(_WIN32)
    #if defined(BALLISTIC_C_EXPORTS)
        #define BALLISTIC_C_API __declspec(dllexport)
    #else
        #define BALLISTIC_C_API __declspec(dllimport)
    #endif
#else
    #define BALLISTIC_C_API __attribute__((visibility("default")))
#endif

#ifdef __cplusplus
extern "C" {
#endif

typedef struct BallisticInputs
{
    double relPos0[3];
    double relVel[3];
    double v0;
    double kDrag;

    double dt;
    double tMax;
    double tolMiss;
    int32_t maxIter;
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

BALLISTIC_C_API int32_t ballistic_solve(const BallisticInputs* in, BallisticOutputs* out);

#ifdef __cplusplus
}
#endif