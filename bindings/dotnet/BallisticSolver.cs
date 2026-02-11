using System;
using System.Numerics;
using System.Runtime.InteropServices;
using System.Text;

public static class Ballistic
{
    public struct Input
    {
        public Vector3 RelPos0;
        public Vector3 RelVel;
        public Vector3 Wind;
        public double V0;
        public double KDrag;

        public int ArcMode;   // 0=Low, 1=High
        public double G;
        public double Dt;
        public double TMax;
        public double TolMiss;
        public int MaxIter;
    }

    public readonly struct Output
    {
        public readonly int CallRc;
        public readonly bool Success;
        public readonly int Status;
        public readonly double Theta, Phi, Miss, TStar;
        public readonly Vector3 RelMissAtStar;
        public readonly string Message;

        internal Output(int rc, in NativeOutputs no, string msg)
        {
            CallRc = rc;
            Success = no.success != 0;
            Status = no.status;
            Theta = no.theta;
            Phi = no.phi;
            Miss = no.miss;
            TStar = no.tStar;
            RelMissAtStar = new Vector3((float)no.relMissAtStar[0], (float)no.relMissAtStar[1], (float)no.relMissAtStar[2]);
            Message = msg;
        }
    }

    public readonly struct ClosestApproachResult
    {
        public readonly int CallRc;
        public readonly double TStar;
        public readonly double Miss;
        public readonly Vector3 RelMissAtStar;

        internal ClosestApproachResult(int rc, double tStar, double miss, Vector3 relMiss)
        {
            CallRc = rc;
            TStar = tStar;
            Miss = miss;
            RelMissAtStar = relMiss;
        }
    }

    public readonly struct VacuumArcResult
    {
        public readonly int CallRc;     // <0: API error
        public readonly bool Reachable; // true if returns 1
        public readonly double Theta;
        public readonly double Phi;

        internal VacuumArcResult(int rc, bool reachable, double theta, double phi)
        {
            CallRc = rc;
            Reachable = reachable;
            Theta = theta;
            Phi = phi;
        }
    }

    public static Input InputDefaults()
    {
        NativeInputs ni = default;
        Native.ballistic_inputs_init(ref ni);

        return new Input
        {
            RelPos0 = Vector3.Zero,
            RelVel = Vector3.Zero,
            Wind = new Vector3((float)ni.wind[0], (float)ni.wind[1], (float)ni.wind[2]),
            V0 = 0.0,
            KDrag = 0.0,
            ArcMode = ni.arcMode,
            G = ni.g,
            Dt = ni.dt,
            TMax = ni.tMax,
            TolMiss = ni.tolMiss,
            MaxIter = ni.maxIter,
        };
    }

    public static Output Solve(Input input)
    {
        NativeInputs ni = default;
        Native.ballistic_inputs_init(ref ni);

        FillNativeInputs(ref ni, input);

        int rc = Native.ballistic_solve(ref ni, out NativeOutputs no);

        string msg;
        unsafe
        {
            fixed (byte* p = no.message) msg = CStr(p, 256);
        }

        return new Output(rc, in no, msg);
    }

    public static void Rk4Step(ref Vector3 r, ref Vector3 v, double h, double kDrag, double g, Vector3 wind)
    {
        double[] r3 = { r.X, r.Y, r.Z };
        double[] v3 = { v.X, v.Y, v.Z };
        double[] w3 = { wind.X, wind.Y, wind.Z };

        Native.ballistic_rk4_step(r3, v3, h, kDrag, g, w3);

        r = new Vector3((float)r3[0], (float)r3[1], (float)r3[2]);
        v = new Vector3((float)v3[0], (float)v3[1], (float)v3[2]);
    }

    public static (double[] t, Vector3[] r, Vector3[] v) SimulateTrajectory(
        Vector3 r0,
        Vector3 v0,
        double kDrag,
        double g,
        Vector3 wind,
        double dt,
        int steps)
    {
        if (steps < 0) throw new ArgumentOutOfRangeException(nameof(steps));

        int count = steps + 1;

        double[] pos3 = new double[3 * count];
        double[] vel3 = new double[3 * count];
        double[] t = new double[count];

        double[] r0_3 = { r0.X, r0.Y, r0.Z };
        double[] v0_3 = { v0.X, v0.Y, v0.Z };
        double[] w3 = { wind.X, wind.Y, wind.Z };

        int outCount = 0;
        int rc = Native.ballistic_simulate_trajectory(
            r0_3, v0_3,
            kDrag, g, w3,
            dt, steps,
            pos3, vel3, t, ref outCount);

        if (rc != 0)
            throw new InvalidOperationException($"ballistic_simulate_trajectory failed: rc={rc}");

        Vector3[] r = new Vector3[outCount];
        Vector3[] v = new Vector3[outCount];

        for (int i = 0; i < outCount; ++i)
        {
            int j = 3 * i;
            r[i] = new Vector3((float)pos3[j + 0], (float)pos3[j + 1], (float)pos3[j + 2]);
            v[i] = new Vector3((float)vel3[j + 0], (float)vel3[j + 1], (float)vel3[j + 2]);
        }

        if (outCount != count)
        {
            Array.Resize(ref t, outCount);
        }

        return (t, r, v);
    }

    public static (double[] t, Vector3[] r, Vector3[] v) SimulateTrajectoryFromAngles(
        Vector3 r0,
        double speed,
        double theta,
        double phi,
        double kDrag,
        double g,
        Vector3 wind,
        double dt,
        int steps)
    {
        if (steps < 0) throw new ArgumentOutOfRangeException(nameof(steps));

        int count = steps + 1;

        double[] pos3 = new double[3 * count];
        double[] vel3 = new double[3 * count];
        double[] t = new double[count];

        double[] r0_3 = { r0.X, r0.Y, r0.Z };
        double[] w3 = { wind.X, wind.Y, wind.Z };

        int outCount = 0;
        int rc = Native.ballistic_simulate_trajectory_from_angles(
            r0_3,
            speed, theta, phi,
            kDrag, g, w3,
            dt, steps,
            pos3, vel3, t, ref outCount);

        if (rc != 0)
            throw new InvalidOperationException($"ballistic_simulate_trajectory_from_angles failed: rc={rc}");

        Vector3[] r = new Vector3[outCount];
        Vector3[] v = new Vector3[outCount];

        for (int i = 0; i < outCount; ++i)
        {
            int j = 3 * i;
            r[i] = new Vector3((float)pos3[j + 0], (float)pos3[j + 1], (float)pos3[j + 2]);
            v[i] = new Vector3((float)vel3[j + 0], (float)vel3[j + 1], (float)vel3[j + 2]);
        }

        if (outCount != count)
        {
            Array.Resize(ref t, outCount);
        }

        return (t, r, v);
    }

    public static ClosestApproachResult ClosestApproach(Input input, double theta, double phi)
    {
        NativeInputs ni = default;
        Native.ballistic_inputs_init(ref ni);
        FillNativeInputs(ref ni, input);

        double tStar = 0.0;
        double miss = 0.0;
        double[] relMiss3 = new double[3];

        int rc = Native.ballistic_find_closest_approach(ref ni, theta, phi, ref tStar, relMiss3, ref miss);

        Vector3 relMiss = new Vector3((float)relMiss3[0], (float)relMiss3[1], (float)relMiss3[2]);
        return new ClosestApproachResult(rc, tStar, miss, relMiss);
    }

    public static VacuumArcResult VacuumArcAnglesToPoint(Vector3 R, double v0, int arcMode, double g)
    {
        double[] R3 = { R.X, R.Y, R.Z };
        double theta = 0.0, phi = 0.0;

        int rc = Native.ballistic_vacuum_arc_angles_to_point(R3, v0, arcMode, g, ref theta, ref phi);

        if (rc < 0) return new VacuumArcResult(rc, false, theta, phi);
        return new VacuumArcResult(0, rc == 1, theta, phi);
    }

    public static (double theta, double phi) VacuumLeadInitialGuess(Vector3 relPos0, Vector3 relVel, double v0, int arcMode, double g)
    {
        double[] rp = { relPos0.X, relPos0.Y, relPos0.Z };
        double[] rv = { relVel.X, relVel.Y, relVel.Z };

        double theta = 0.0, phi = 0.0;
        Native.ballistic_initial_guess_vacuum_lead(rp, rv, v0, arcMode, g, ref theta, ref phi);
        return (theta, phi);
    }

    private static void FillNativeInputs(ref NativeInputs ni, in Input input)
    {
        unsafe
        {
            ni.relPos0[0] = input.RelPos0.X; ni.relPos0[1] = input.RelPos0.Y; ni.relPos0[2] = input.RelPos0.Z;
            ni.relVel[0] = input.RelVel.X; ni.relVel[1] = input.RelVel.Y; ni.relVel[2] = input.RelVel.Z;
            ni.wind[0] = input.Wind.X; ni.wind[1] = input.Wind.Y; ni.wind[2] = input.Wind.Z;
        }

        ni.v0 = input.V0;
        ni.kDrag = input.KDrag;

        ni.arcMode = input.ArcMode;
        ni.g = input.G;
        ni.dt = input.Dt;
        ni.tMax = input.TMax;
        ni.tolMiss = input.TolMiss;
        ni.maxIter = input.MaxIter;
    }

    private static unsafe string CStr(byte* p, int maxBytes)
    {
        int n = 0;
        for (; n < maxBytes; n++) if (p[n] == 0) break;
        return Encoding.UTF8.GetString(p, n);
    }

    // ---------- Native layer ----------
    [StructLayout(LayoutKind.Sequential, Pack = 8)]
    private unsafe struct NativeInputs
    {
        public fixed double relPos0[3];
        public fixed double relVel[3];
        public double v0;
        public double kDrag;
        public int arcMode;
        public int _pad0;
        public double g;
        public fixed double wind[3];
        public double dt;
        public double tMax;
        public double tolMiss;
        public int maxIter;
        public int _pad1;
    }

    [StructLayout(LayoutKind.Sequential)]
    private unsafe struct NativeOutputs
    {
        public int success;
        public int status;
        public double theta;
        public double phi;
        public double miss;
        public double tStar;
        public fixed double relMissAtStar[3];
        public fixed byte message[256];
    }

    private static class Native
    {
        private const string Dll = "ballistic_solver";

        [DllImport(Dll, CallingConvention = CallingConvention.Cdecl, EntryPoint = "ballistic_inputs_init")]
        public static extern void ballistic_inputs_init(ref NativeInputs input);

        [DllImport(Dll, CallingConvention = CallingConvention.Cdecl, EntryPoint = "ballistic_solve")]
        public static extern int ballistic_solve(ref NativeInputs input, out NativeOutputs output);

        [DllImport(Dll, CallingConvention = CallingConvention.Cdecl, EntryPoint = "ballistic_rk4_step")]
        public static extern void ballistic_rk4_step(
            [In, Out] double[] r3,
            [In, Out] double[] v3,
            double h,
            double kDrag,
            double g,
            [In] double[] wind3);

        [DllImport(Dll, CallingConvention = CallingConvention.Cdecl, EntryPoint = "ballistic_simulate_trajectory")]
        public static extern int ballistic_simulate_trajectory(
            [In] double[] r0_3,
            [In] double[] v0_3,
            double kDrag,
            double g,
            [In] double[] wind3,
            double dt,
            int steps,
            [Out] double[] out_pos3,
            [Out] double[] out_vel3,
            [Out] double[] out_t,
            ref int out_count);

        [DllImport(Dll, CallingConvention = CallingConvention.Cdecl, EntryPoint = "ballistic_simulate_trajectory_from_angles")]
        public static extern int ballistic_simulate_trajectory_from_angles(
            [In] double[] r0_3,
            double speed,
            double theta,
            double phi,
            double kDrag,
            double g,
            [In] double[] wind3,
            double dt,
            int steps,
            [Out] double[] out_pos3,
            [Out] double[] out_vel3,
            [Out] double[] out_t,
            ref int out_count);

        [DllImport(Dll, CallingConvention = CallingConvention.Cdecl, EntryPoint = "ballistic_find_closest_approach")]
        public static extern int ballistic_find_closest_approach(
            ref NativeInputs input,
            double theta,
            double phi,
            ref double out_tStar,
            [Out] double[] out_relMissAtStar3,
            ref double out_miss);

        [DllImport(Dll, CallingConvention = CallingConvention.Cdecl, EntryPoint = "ballistic_vacuum_arc_angles_to_point")]
        public static extern int ballistic_vacuum_arc_angles_to_point(
            [In] double[] R3,
            double v0,
            int arcMode,
            double g,
            ref double out_theta,
            ref double out_phi);

        [DllImport(Dll, CallingConvention = CallingConvention.Cdecl, EntryPoint = "ballistic_initial_guess_vacuum_lead")]
        public static extern void ballistic_initial_guess_vacuum_lead(
            [In] double[] relPos0_3,
            [In] double[] relVel_3,
            double v0,
            int arcMode,
            double g,
            ref double out_theta,
            ref double out_phi);
    }
}
