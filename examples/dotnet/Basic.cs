using System;
using System.Numerics;

class Program
{
    static void Main()
    {
        var input = Ballistic.InputDefaults();
        input.RelPos0 = new Vector3(100, 30, 10);
        input.RelVel  = new Vector3(-10, 30, 0);
        input.V0 = 80;
        input.KDrag = 0.005;

        var output = Ballistic.Solve(input);

        Console.WriteLine($"call_ok={output.CallRc}");
        Console.WriteLine($"success={output.Success}");
        Console.WriteLine($"status={output.Status}");
        Console.WriteLine($"theta={output.Theta}");
        Console.WriteLine($"phi={output.Phi}");
        Console.WriteLine($"miss={output.Miss}");
        Console.WriteLine($"tStar={output.TStar}");
        Console.WriteLine($"relMiss=[{output.RelMissAtStar.X}, {output.RelMissAtStar.Y}, {output.RelMissAtStar.Z}]");
        Console.WriteLine($"message={output.Message}");
    }
}
