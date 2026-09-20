// Minimal use of the public C++ surface. See <ballistic_solver.hpp>.

#include <chrono>
#include <iostream>

#include "ballistic_solver.hpp"

int main()
{
    // Target state relative to the projectile at t = 0
    bs::Problem problem;
    problem.rel_pos0 = { 100.0, 30.0, 10.0 };
    problem.rel_vel  = { -10.0, 30.0, 0.0 };

    problem.v0     = 80.0;     // muzzle speed
    problem.k_drag = 0.005;    // quadratic drag coefficient

    // Solver configuration. Every field has a usable default; the numeric
    // overrides stay at 0 to mean "use the preset's value".
    bs::Options options;
    //options.preset   = bs::Preset::Precise;
    //options.arc      = bs::Arc::High;
    //options.wind     = { 0.0, 0.0, 0.0 };
    //options.gravity  = 9.80665;
    //options.tol_miss = 1e-2;

    const auto t0 = std::chrono::steady_clock::now();
    const bs::Intercept r = bs::solve(problem, options);
    const auto t1 = std::chrono::steady_clock::now();

    constexpr double rad_to_deg = 180.0 / 3.141592653589793238462643383279502884;

    std::cout << "success    : " << r.success << "\n";
    std::cout << "elevation  : " << r.theta * rad_to_deg << " deg\n";
    std::cout << "azimuth    : " << r.phi * rad_to_deg << " deg\n";
    std::cout << "miss       : " << r.miss << "\n";
    std::cout << "time       : " << r.t_star << "\n";
    std::cout << "iterations : " << r.iterations << "\n";
    std::cout << "message    : " << r.message
              << " (status : " << static_cast<int>(r.status) << ")\n";

    const double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
    std::cout << "elapsed    : " << ms << " ms\n";

    return r ? 0 : 1;
}
