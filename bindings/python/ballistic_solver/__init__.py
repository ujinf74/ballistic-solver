from ._core import (
    ArcMode,
    BallisticParams,
    closest_approach_accel,
    k_drag_from_physical,
    params_preset,
    solve,
    solve_accel,
    rk4_step,
    simulate,
    closest_approach,
    vacuum_arc_angles_to_point,
    vacuum_lead_initial_guess,
)
from ._version import version as __version__


class Solver:
    def __init__(self, params=None, arcMode=None):
        self.params = params if params is not None else BallisticParams()
        self.arcMode = arcMode

    @classmethod
    def preset(cls, preset="balanced", arcMode=None):
        return cls(params_preset(preset), arcMode=arcMode)

    def solve(self, relPos0, relVel, v0, kDrag, relAcc=None):
        return solve(relPos0, relVel, v0, kDrag, self.arcMode, self.params, relAcc)

    def solve_accel(self, relPos0, relVel, relAcc, v0, kDrag):
        return solve_accel(relPos0, relVel, relAcc, v0, kDrag, self.arcMode, self.params)


class Projectile:
    def __init__(self, mass, area, dragCoefficient):
        self.mass = mass
        self.area = area
        self.dragCoefficient = dragCoefficient

    @property
    def kDrag(self):
        return k_drag_from_physical(1.225, self.dragCoefficient, self.area, self.mass)

    def k_drag(self, airDensity=1.225):
        return k_drag_from_physical(airDensity, self.dragCoefficient, self.area, self.mass)


class Environment:
    def __init__(self, airDensity=1.225, wind=(0.0, 0.0, 0.0), g=9.80665):
        self.airDensity = airDensity
        self.wind = wind
        self.g = g

    def params(self, preset="balanced"):
        p = params_preset(preset)
        p.wind = self.wind
        p.g = self.g
        return p


__all__ = [
    "ArcMode",
    "BallisticParams",
    "Environment",
    "Projectile",
    "Solver",
    "closest_approach_accel",
    "k_drag_from_physical",
    "params_preset",
    "solve",
    "solve_accel",
    "rk4_step",
    "simulate",
    "closest_approach",
    "vacuum_arc_angles_to_point",
    "vacuum_lead_initial_guess",
]
