""" 87Rb

References:
Daniel A. Steck, ‘‘Rubidium 87 D Line Data,’’ available online at http://steck.us/alkalidata (revision 2.3.2, 10 September 2023). """

import numpy as np
import typing
import atomic_physics as ap
import scipy.constants as consts


# level aliases
ground_level = S12 = ap.Level(n=5, S=1 / 2, L=0, J=1 / 2)
P12 = ap.Level(n=5, S=1 / 2, L=1, J=1 / 2)
P32 = ap.Level(n=5, S=1 / 2, L=1, J=3 / 2)


class Rb87(ap.Atom):
    def __init__(
        self,
        *,
        B: typing.Optional[float] = None,
        level_filter: typing.Optional[typing.List[ap.Level]] = None
    ):
        """Rb87 atomic structure.

        :param B: B-field in Tesla (can be changed using :meth setB:)
        :param level_filter: list of Levels to include in the simulation, if
            None we include all levels.
        """
        levels = {
            ground_level: ap.LevelData(
                Ahfs=3.417341305452155e9 * consts.h, 
                g_J=2.0023311320, 
                g_I=-0.000955141410, 
            ),
            P12: ap.LevelData(
                Ahfs=408.32815e6 * consts.h,
                g_I=-0.000955141410,
                g_J=0.666,
            ),
            P32: ap.LevelData(
                Ahfs=84.4717530e6 * consts.h,
                Bhfs=12.496537e6*consts.h,
                g_I=-0.000955141410,
                g_J=1.336213,
            ),

        }

        transitions = {
            "780": ap.Transition(
                lower=ground_level,
                upper=P12,
                A=1/(26.244*10**(-9)),
                freq=2 * np.pi * 384.230484468562*10**9, 
            ),
            "795": ap.Transition(
                lower=ground_level,
                upper=P32,
                A=1/(27.704*10**(-9)),
                freq=2 * np.pi * 377.170746354*10**9, 
            ),
        }

        super().__init__(
            B=B, I=3/2, levels=levels, transitions=transitions, level_filter=level_filter
        )