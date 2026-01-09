using System;
using System.Runtime.InteropServices;

[StructLayout(LayoutKind.Sequential, Pack = 8)]
public struct BallisticInputs
{
    [MarshalAs(UnmanagedType.ByValArray, SizeConst = 3)]
    public double[] relPos0;

    [MarshalAs(UnmanagedType.ByValArray, SizeConst = 3)]
    public double[] relVel;

    public double v0;
    public double kDrag;

    public double dt;
    public double tMax;
    public double tolMiss;
    public Int32 maxIter;
}

[StructLayout(LayoutKind.Sequential, Pack = 8)]
public struct BallisticOutputs
{
    public Int32 success;
    public Int32 status;

    public double theta;
    public double phi;
    public double miss;
    public double tStar;

    [MarshalAs(UnmanagedType.ByValArray, SizeConst = 3)]
    public double[] relMissAtStar;

    [MarshalAs(UnmanagedType.ByValTStr, SizeConst = 256)]
    public string message;
}

public static class BallisticNative
{
    // Windows: ballistic_c.dll must be in the executable directory or PATH
    [DllImport("ballistic_c.dll", CallingConvention = CallingConvention.Cdecl)]
    public static extern Int32 ballistic_solve(ref BallisticInputs input, out BallisticOutputs output);
}

public class Program
{
    public static void Main()
    {
        var input = new BallisticInputs
        {
            relPos0 = new double[] { 120.0, 30.0, 5.0 },
            relVel  = new double[] { 2.0, -1.0, 0.0 },
            v0 = 90.0,
            kDrag = 0.002,
            dt = 0.01,
            tMax = 30.0,
            tolMiss = 0.5,
            maxIter = 30
        };

        BallisticOutputs output;
        int ok = BallisticNative.ballistic_solve(ref input, out output);

        Console.WriteLine($"call_ok={ok}");
        Console.WriteLine($"success={output.success}");
        Console.WriteLine($"status={output.status}");
        Console.WriteLine($"theta={output.theta}");
        Console.WriteLine($"phi={output.phi}");
        Console.WriteLine($"miss={output.miss}");
        Console.WriteLine($"tStar={output.tStar}");
        Console.WriteLine($"relMiss=[{output.relMissAtStar[0]}, {output.relMissAtStar[1]}, {output.relMissAtStar[2]}]");
        Console.WriteLine($"message={output.message}");
    }
}