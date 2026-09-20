/* Smoke test for the C ABI declared in <ballistic_solver_c_api.h>.

   Uses an explicit CHECK instead of assert(): CI builds Release, which defines
   NDEBUG and would compile every assert away -- including any API call written
   inside one. Calls are made on their own line and only their result is
   checked. Failures are counted and reported through the exit code. */

#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>

#include "ballistic_solver_c_api.h"

static int g_failures = 0;

static void check(int ok, const char* expr, int line)
{
    if (!ok)
    {
        fprintf(stderr, "FAIL (line %d): %s\n", line, expr);
        ++g_failures;
    }
}

#define CHECK(expr) check((expr) ? 1 : 0, #expr, __LINE__)

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

    CHECK(ballistic_solver_abi_version() == BALLISTIC_SOLVER_ABI_VERSION);

    {
        int32_t rc = ballistic_solve(&in, &out);
        CHECK(rc == 0);
        CHECK(out.success != 0);
        CHECK(out.status == 0);
    }

    CHECK(isfinite(out.theta));
    CHECK(isfinite(out.phi));
    CHECK(isfinite(out.miss));
    CHECK(isfinite(out.tStar));
    CHECK(is_finite3(out.relMissAtStar));
    CHECK(out.miss >= 0.0);
    CHECK(out.iterations > 0);
    CHECK(out.acceptedSteps >= 0);
    /* lastLambda/lastAlpha are Levenberg-Marquardt diagnostics. The default
       core (coordinate-residual) uses neither damping nor a line search, so it
       leaves them at NaN -- do not assert they are finite here. */

    printf("abi=%u\n", ballistic_solver_abi_version());
    printf("version=%s\n", ballistic_solver_version_string());

    printf("success=%d\n", out.success);
    printf("status=%d\n", out.status);
    printf("theta=%.17g rad\n", out.theta);
    printf("phi=%.17g rad\n", out.phi);
    printf("miss=%.17g m\n", out.miss);
    printf("t*=%.17g s\n", out.tStar);
    printf("iterations=%d\n", out.iterations);
    printf("acceptedSteps=%d\n", out.acceptedSteps);
    printf("lastLambda=%.17g\n", out.lastLambda);
    printf("lastAlpha=%.17g\n", out.lastAlpha);
    printf("relMissAtStar=[%.17g, %.17g, %.17g]\n",
        out.relMissAtStar[0], out.relMissAtStar[1], out.relMissAtStar[2]);
    printf("message=%s\n", out.message);

    {
        double k = 0.0;
        int32_t rc = ballistic_k_drag_from_physical(1.225, 0.30, 0.00426, 0.145, &k);
        CHECK(rc == 0);
        CHECK(isfinite(k));
        CHECK(k > 0.0);
    }

    {
        BallisticInputs precise = in;
        int32_t rc = ballistic_inputs_apply_preset(&precise, 2);
        CHECK(rc == 0);
        CHECK(precise.tolMiss < in.tolMiss);
        CHECK(precise.maxIter > in.maxIter);
        CHECK(precise.preset == 2); /* preset must be recorded so solve re-expands full tuning */

        /* Solve with the precise preset: deep tuning (line-search/lambda/FD/golden)
           now flows through ballistic_solve, so the tighter tolerance is reachable. */
        BallisticOutputs preciseOut;
        memset(&preciseOut, 0, sizeof(preciseOut));

        rc = ballistic_solve(&precise, &preciseOut);
        CHECK(rc == 0);
        CHECK(preciseOut.success != 0);
        CHECK(preciseOut.miss <= precise.tolMiss);

        printf("precise: miss=%.17g m (tol=%.17g)\n", preciseOut.miss, precise.tolMiss);
    }

    {
        BallisticAccelInputs accIn;
        BallisticOutputs accOut;
        int32_t rc;

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

        rc = ballistic_solve_accel(&accIn, &accOut);
        CHECK(rc == 0);
        CHECK(accOut.success != 0);
        CHECK(isfinite(accOut.theta));
        CHECK(isfinite(accOut.phi));
        CHECK(isfinite(accOut.miss));

        printf("accel: miss=%.17g m t*=%.17g s\n", accOut.miss, accOut.tStar);
    }

    if (g_failures != 0)
    {
        fprintf(stderr, "smoke_test: %d check(s) failed\n", g_failures);
        return 1;
    }

    printf("smoke_test: ok\n");
    return 0;
}
