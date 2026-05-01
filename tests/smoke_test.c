#include <assert.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>

#include "ballistic_solver_c_api.h"

static int is_finite3(const double v[3])
{
    return isfinite(v[0]) && isfinite(v[1]) && isfinite(v[2]);
}

int main(void)
{
    BallisticInputs in;
    BallisticOutputs out;

    ballistic_inputs_init(&in);
    memset(&out, 0, sizeof(out));

    in.relPos0[0] = 120.0;
    in.relPos0[1] = 30.0;
    in.relPos0[2] = 5.0;

    in.relVel[0] = 2.0;
    in.relVel[1] = -1.0;
    in.relVel[2] = 0.0;

    in.v0 = 90.0;
    in.kDrag = 0.002;

    assert(ballistic_solver_abi_version() == BALLISTIC_SOLVER_ABI_VERSION);

    {
        int32_t rc = ballistic_solve(&in, &out);
        assert(rc == 0);
    }

    assert(isfinite(out.theta));
    assert(isfinite(out.phi));
    assert(isfinite(out.miss));
    assert(isfinite(out.tStar));
    assert(is_finite3(out.relMissAtStar));
    assert(out.miss >= 0.0);

    printf("abi=%u\n", ballistic_solver_abi_version());
    printf("version=%s\n", ballistic_solver_version_string());

    printf("success=%d\n", out.success);
    printf("status=%d\n", out.status);
    printf("theta=%.17g rad\n", out.theta);
    printf("phi=%.17g rad\n", out.phi);
    printf("miss=%.17g m\n", out.miss);
    printf("t*=%.17g s\n", out.tStar);
    printf("relMissAtStar=[%.17g, %.17g, %.17g]\n",
        out.relMissAtStar[0], out.relMissAtStar[1], out.relMissAtStar[2]);
    printf("message=%s\n", out.message);

    {
        double k = 0.0;
        assert(ballistic_k_drag_from_physical(1.225, 0.30, 0.00426, 0.145, &k) == 0);
        assert(isfinite(k));
        assert(k > 0.0);
    }

    {
        BallisticInputs precise;
        ballistic_inputs_init(&precise);
        assert(ballistic_inputs_apply_preset(&precise, 2) == 0);
        assert(precise.tolMiss < in.tolMiss);
        assert(precise.maxIter > in.maxIter);
    }

    {
        BallisticAccelInputs accIn;
        BallisticOutputs accOut;
        ballistic_accel_inputs_init(&accIn);
        memset(&accOut, 0, sizeof(accOut));

        accIn.base.relPos0[0] = 120.0;
        accIn.base.relPos0[1] = 30.0;
        accIn.base.relPos0[2] = 5.0;
        accIn.base.relVel[0] = 2.0;
        accIn.base.relVel[1] = -1.0;
        accIn.base.relVel[2] = 0.0;
        accIn.relAcc[0] = 0.0;
        accIn.relAcc[1] = 0.2;
        accIn.relAcc[2] = 0.0;
        accIn.base.v0 = 90.0;
        accIn.base.kDrag = 0.002;

        assert(ballistic_solve_accel(&accIn, &accOut) == 0);
        assert(isfinite(accOut.theta));
        assert(isfinite(accOut.phi));
        assert(isfinite(accOut.miss));
    }

    return 0;
}
