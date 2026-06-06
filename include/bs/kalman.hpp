#pragma once

#include "predictor.hpp"

// ================================================================
// Constant-acceleration Kalman tracker (position-only input)
// ================================================================
// A position-only target tracker: it ingests (timestamp, relative position)
// measurements and maintains a constant-acceleration state estimate
// (position, velocity, acceleration) with covariance. It produces a
// TargetPredictor for lead-fire, so the caller never supplies relVel/relAcc.
//
// Each Cartesian axis is an independent 1D constant-acceleration filter
// (state [p, v, a]) with a white-noise-jerk process model. This denoises the
// track and recovers velocity/acceleration without finite-difference noise;
// the lead prediction is the second-order extrapolation p + v*tau + 0.5*a*tau^2,
// which a curved (maneuvering) target only follows approximately -- see
// solve_launch_angles_predicted for how that approximation is used.

struct CaKalman1
{
    double x[3] = { 0.0, 0.0, 0.0 }; // [p, v, a]
    double P[3][3] = {
        { 1.0, 0.0, 0.0 },
        { 0.0, 1.0, 0.0 },
        { 0.0, 0.0, 1.0 },
    };

    void reset(double p0, double posVar, double velVar, double accVar)
    {
        x[0] = p0;
        x[1] = 0.0;
        x[2] = 0.0;
        for (int i = 0; i < 3; ++i)
        {
            for (int j = 0; j < 3; ++j)
            {
                P[i][j] = 0.0;
            }
        }
        P[0][0] = posVar;
        P[1][1] = velVar;
        P[2][2] = accVar;
    }

    // Time update over dt with white-noise-jerk spectral density q.
    void predict(double dt, double q)
    {
        const double d2 = 0.5 * dt * dt;

        // x = F x, F = [[1,dt,dt^2/2],[0,1,dt],[0,0,1]]
        const double p = x[0] + dt * x[1] + d2 * x[2];
        const double v = x[1] + dt * x[2];
        const double a = x[2];
        x[0] = p;
        x[1] = v;
        x[2] = a;

        // P = F P F^T
        double FP[3][3];
        for (int j = 0; j < 3; ++j)
        {
            FP[0][j] = P[0][j] + dt * P[1][j] + d2 * P[2][j];
            FP[1][j] = P[1][j] + dt * P[2][j];
            FP[2][j] = P[2][j];
        }
        double Pn[3][3];
        for (int i = 0; i < 3; ++i)
        {
            Pn[i][0] = FP[i][0] + dt * FP[i][1] + d2 * FP[i][2];
            Pn[i][1] = FP[i][1] + dt * FP[i][2];
            Pn[i][2] = FP[i][2];
        }

        // Q (white-noise-jerk), spectral density q
        const double dt2 = dt * dt;
        const double dt3 = dt2 * dt;
        const double dt4 = dt3 * dt;
        const double dt5 = dt4 * dt;
        const double Q[3][3] = {
            { q * dt5 / 20.0, q * dt4 / 8.0, q * dt3 / 6.0 },
            { q * dt4 / 8.0,  q * dt3 / 3.0, q * dt2 / 2.0 },
            { q * dt3 / 6.0,  q * dt2 / 2.0, q * dt },
        };

        for (int i = 0; i < 3; ++i)
        {
            for (int j = 0; j < 3; ++j)
            {
                P[i][j] = Pn[i][j] + Q[i][j];
            }
        }
    }

    // Measurement update with position z and measurement variance r. H = [1,0,0].
    void update(double z, double r)
    {
        const double S = P[0][0] + r;
        if (!(S > 0.0))
        {
            return;
        }
        const double K0 = P[0][0] / S;
        const double K1 = P[1][0] / S;
        const double K2 = P[2][0] / S;

        const double y = z - x[0];
        x[0] += K0 * y;
        x[1] += K1 * y;
        x[2] += K2 * y;

        // P = (I - K H) P  ->  P[i][j] -= K[i] * P[0][j]
        const double row0[3] = { P[0][0], P[0][1], P[0][2] };
        const double K[3] = { K0, K1, K2 };
        for (int i = 0; i < 3; ++i)
        {
            for (int j = 0; j < 3; ++j)
            {
                P[i][j] -= K[i] * row0[j];
            }
        }
    }
};

struct TargetTracker
{
    CaKalman1 axis[3];
    double q = 1.0;       // process noise (jerk spectral density)
    double r = 1.0;       // measurement variance (position)
    double lastT = 0.0;
    bool started = false;

    double velVar0 = 1.0e3;
    double accVar0 = 1.0e3;

    explicit TargetTracker(double processNoise = 1.0, double measNoise = 1.0)
        : q(processNoise), r(measNoise)
    {
    }

    void update(double t, const Vec3& z)
    {
        const double zc[3] = { z.x, z.y, z.z };
        if (!started)
        {
            for (int i = 0; i < 3; ++i)
            {
                axis[i].reset(zc[i], r, velVar0, accVar0);
            }
            lastT = t;
            started = true;
            return;
        }

        const double dt = t - lastT;
        for (int i = 0; i < 3; ++i)
        {
            if (dt > 0.0)
            {
                axis[i].predict(dt, q);
            }
            axis[i].update(zc[i], r);
        }
        lastT = t;
    }

    Vec3 position() const { return { axis[0].x[0], axis[1].x[0], axis[2].x[0] }; }
    Vec3 velocity() const { return { axis[0].x[1], axis[1].x[1], axis[2].x[1] }; }
    Vec3 acceleration() const { return { axis[0].x[2], axis[1].x[2], axis[2].x[2] }; }

    // Predicted relative position tau seconds ahead of the last measurement.
    Vec3 predict(double tau) const
    {
        Vec3 out{};
        out.x = axis[0].x[0] + tau * axis[0].x[1] + 0.5 * tau * tau * axis[0].x[2];
        out.y = axis[1].x[0] + tau * axis[1].x[1] + 0.5 * tau * tau * axis[1].x[2];
        out.z = axis[2].x[0] + tau * axis[2].x[1] + 0.5 * tau * tau * axis[2].x[2];
        return out;
    }

    // A TargetPredictor for solve_launch_angles_predicted (snapshot of current state).
    TargetPredictor predictor() const
    {
        const Vec3 p = position();
        const Vec3 v = velocity();
        const Vec3 a = acceleration();
        TargetPredictor tp;
        tp.pos = [p, v, a](double tau) -> Vec3
        {
            return { p.x + tau * v.x + 0.5 * tau * tau * a.x,
                     p.y + tau * v.y + 0.5 * tau * tau * a.y,
                     p.z + tau * v.z + 0.5 * tau * tau * a.z };
        };
        return tp;
    }
};
