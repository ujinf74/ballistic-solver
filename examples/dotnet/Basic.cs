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

        Console.WriteLine($"call_rc={output.CallRc}");
        Console.WriteLine($"success={output.Success}");
        Console.WriteLine($"status={output.Status}");
        if (output.CallRc != 0 || !output.Success)
        {
            Console.WriteLine($"solve failed: rc={output.CallRc} status={output.Status} message={output.Message}");
            Environment.ExitCode = 1;
            return;
        }
        Console.WriteLine($"theta={output.Theta}");
        Console.WriteLine($"phi={output.Phi}");
        Console.WriteLine($"miss={output.Miss}");
        Console.WriteLine($"tStar={output.TStar}");
        Console.WriteLine($"iterations={output.Iterations}");
        Console.WriteLine($"acceptedSteps={output.AcceptedSteps}");
        Console.WriteLine($"lastLambda={output.LastLambda}");
        Console.WriteLine($"lastAlpha={output.LastAlpha}");
        Console.WriteLine($"relMiss=[{output.RelMissAtStar.X}, {output.RelMissAtStar.Y}, {output.RelMissAtStar.Z}]");
        Console.WriteLine($"message={output.Message}");
    }
}
