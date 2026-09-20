// Smoke test for the public C++ surface declared in <ballistic_solver.hpp>.
//
// Uses an explicit CHECK instead of assert(): CI builds Release, which defines
// NDEBUG and would compile every assert away. Failures are counted and reported
// through the exit code.

#include <cmath>
#include <cstdio>

#include "ballistic_solver.hpp"

namespace {

int g_failures = 0;

void check(bool ok, const char* expr, int line)
{
    if (!ok)
    {
        std::fprintf(stderr, "FAIL (line %d): %s\n", line, expr);
        ++g_failures;
    }
}

#define CHECK(expr) check((expr), #expr, __LINE__)

bool finite3(const bs::Vec3& v)
{
    return std::isfinite(v.x) && std::isfinite(v.y) && std::isfinite(v.z);
}

bs::Problem base_problem()
{
    bs::Problem p;
    p.rel_pos0 = { 120.0, 30.0, 5.0 };
    p.rel_vel  = { 2.0, -1.0, 0.0 };
    p.v0       = 90.0;
    p.k_drag   = 0.002;
    return p;
}

void report(const char* label, const bs::Intercept& r)
{
    std::printf("%-14s success=%d status=%d theta=%.6g phi=%.6g miss=%.6g t*=%.6g iters=%d msg=%.*s\n",
        label, static_cast<int>(r.success), static_cast<int>(r.status),
        r.theta, r.phi, r.miss, r.t_star, r.iterations,
        static_cast<int>(r.message.size()), r.message.data());
}

} // namespace

int main()
{
    // Default solve: coordinate-residual core, Balanced preset (tol_miss 1e-2).
    double lowArcTheta = 0.0;
    {
        const bs::Intercept r = bs::solve(base_problem());
        report("solve", r);

        CHECK(r.success);
        CHECK(static_cast<bool>(r));          // operator bool tracks success
        CHECK(r.status == bs::Status::Ok);
        CHECK(std::isfinite(r.theta));
        CHECK(std::isfinite(r.phi));
        CHECK(std::isfinite(r.miss));
        CHECK(std::isfinite(r.t_star));
        CHECK(r.miss >= 0.0);
        CHECK(r.miss <= 1e-2);                // Balanced preset tolerance
        CHECK(r.t_star > 0.0);
        CHECK(finite3(r.rel_miss_at_star));
        CHECK(r.iterations > 0);
        CHECK(!r.message.empty());

        lowArcTheta = r.theta;
    }

    // Preset must reach the solver: Precise tightens tol_miss to 1e-5.
    {
        bs::Options o;
        o.preset = bs::Preset::Precise;

        const bs::Intercept r = bs::solve(base_problem(), o);
        report("precise", r);

        CHECK(r.success);
        CHECK(r.status == bs::Status::Ok);
        CHECK(r.miss <= 1e-5);
    }

    // Arc must reach the solver: the high arc is a steeper shot than the low one.
    {
        bs::Options o;
        o.arc = bs::Arc::High;

        const bs::Intercept r = bs::solve(base_problem(), o);
        report("high arc", r);

        CHECK(r.success);
        CHECK(r.theta > lowArcTheta);
    }

    // Constant relative acceleration path.
    {
        bs::Problem p = base_problem();
        p.rel_acc = { 0.0, 0.2, 0.0 };

        const bs::Intercept r = bs::solve(p);
        report("accel", r);

        CHECK(r.success);
        CHECK(r.status == bs::Status::Ok);
        CHECK(finite3(r.rel_miss_at_star));
    }

    // Auxiliary-residual variant.
    {
        const bs::Intercept r = bs::solve_aux(base_problem());
        report("solve_aux", r);

        CHECK(r.success);
        CHECK(r.status == bs::Status::Ok);
        CHECK(r.miss <= 1e-2);
    }

    // Rejected input must come back as a failed, diagnosable result.
    {
        bs::Problem p = base_problem();
        p.v0 = 0.0;

        const bs::Intercept r = bs::solve(p);
        report("invalid v0", r);

        CHECK(!r.success);
        CHECK(!static_cast<bool>(r));
        CHECK(r.status == bs::Status::InvalidInput);
        CHECK(!r.message.empty());
    }

    if (g_failures != 0)
    {
        std::fprintf(stderr, "smoke_test_cpp: %d check(s) failed\n", g_failures);
        return 1;
    }

    std::printf("smoke_test_cpp: ok\n");
    return 0;
}
