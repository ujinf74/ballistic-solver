#include <cstdint>
#include <cmath>
#include <cstdio>
#include <cstring>

#include "ballistic/ballistic_solver.hpp"

#if defined(_WIN32)
    #define BALLISTIC_EXPORT __declspec(dllexport)
#else
    #define BALLISTIC_EXPORT __attribute__((visibility("default")))
#endif

extern "C"
{
    struct BallisticInputs
    {
        double relPos0[3];
        double relVel[3];
        double v0;
        double kDrag;

        double dt;
        double tMax;
        double tolMiss;
        int32_t maxIter;
    };

    struct BallisticOutputs
    {
        int32_t success;
        int32_t status;

        double theta;
        double phi;
        double miss;
        double tStar;

        double relMissAtStar[3];

        char message[256];
    };

    BALLISTIC_EXPORT int32_t ballistic_solve(const BallisticInputs* in, BallisticOutputs* out)
    {
        if (in == nullptr || out == nullptr)
        {
            return 0;
        }

        std::memset(out, 0, sizeof(*out));
        std::snprintf(out->message, sizeof(out->message), "Unknown");

        try
        {
            BallisticParams P;

            if (std::isfinite(in->dt) && in->dt > 0.0) { P.dt = in->dt; }
            if (std::isfinite(in->tMax) && in->tMax > 0.0) { P.tMax = in->tMax; }
            if (std::isfinite(in->tolMiss) && in->tolMiss > 0.0) { P.tolMiss = in->tolMiss; }
            if (in->maxIter > 0) { P.maxIter = static_cast<int>(in->maxIter); }

            Vec3 relPos0 = { in->relPos0[0], in->relPos0[1], in->relPos0[2] };
            Vec3 relVel  = { in->relVel[0],  in->relVel[1],  in->relVel[2]  };

            SolverResult r = solve_launch_angles(relPos0, relVel, in->v0, in->kDrag, P);

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

            return 1;
        }
        catch (...)
        {
            out->success = 0;
            out->status = static_cast<int32_t>(SolveStatus::InvalidInput);
            std::snprintf(out->message, sizeof(out->message), "Exception caught in ballistic_solve");
            return 0;
        }
    }
}
