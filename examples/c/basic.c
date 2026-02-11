#include <stdint.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#include "ballistic_solver_c_api.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

int main(void)
{
    BallisticInputs in;
    BallisticOutputs out;

    ballistic_inputs_init(&in);
    memset(&out, 0, sizeof(out));

    in.relPos0[0] = 100.0; in.relPos0[1] = 30.0; in.relPos0[2] = 10.0;
    in.relVel[0]  = -10.0; in.relVel[1]  = 30.0; in.relVel[2]  = 0.0;

    in.v0 = 80.0;
    in.kDrag = 0.005;

    int32_t rc = ballistic_solve(&in, &out);
    if (rc != 0)
    {
        printf("call failed: rc=%d status=%d message=%s\n", (int)rc, (int)out.status, out.message);
        return 1;
    }

    printf("elev=%.6f deg az=%.6f deg miss=%.6f t=%.6f\n",
           out.theta * 180.0 / M_PI, out.phi * 180.0 / M_PI, out.miss, out.tStar);
    return 0;
}
