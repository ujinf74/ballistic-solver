#include <iostream>
#include <chrono>
#include "ballistic_solver_core.hpp"

#ifndef M_PI
#define M_PI 3.141592653589793238462643383279502884
#endif

int main()
{
    // Target state relative to the projectile at t = 0
    Vec3 relPos0 = { 100.0, 30.0, 10.0 };
    Vec3 relVel  = { -10.0, 30.0, 0.0 };

    // Projectile parameters
    double v0 = 80.0;     // initial speed magnitude
    double kDrag = 0.005; // drag coefficient

    // Solver configuration
    BallisticParams P;
    //P.arcMode = ArcMode::Low;
    //P.g = 9.80665;
    //P.wind = { 0.0, 0.0, 0.0 };
    //P.dt = 0.01;
    //P.tMax = 20.0;
    //P.tolMiss = 1e-2;

    // Solve launch angles
    auto t0 = std::chrono::steady_clock::now();
    SolverResult r = solve_launch_angles(relPos0, relVel, v0, kDrag, P);
    auto t1 = std::chrono::steady_clock::now();

    // Output results
    std::cout << "success    : " << r.success << "\n";
    std::cout << "elevation  : " << r.theta * 180.0 / M_PI << " deg\n";
    std::cout << "azimuth    : " << r.phi * 180.0 / M_PI << " deg\n";
    std::cout << "miss       : " << r.miss << "\n";
    std::cout << "time       : " << r.tStar << "\n";
    std::cout << "iterations : " << r.report.iterations
              << " (accepted : " << r.report.acceptedSteps << ")\n";
    std::cout << "message    : " << r.report.message
              << " (status : " << static_cast<int>(r.report.status) << ")\n";

    double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
    std::cout << "elapsed    : " << ms << " ms\n";

    return 0;
}
