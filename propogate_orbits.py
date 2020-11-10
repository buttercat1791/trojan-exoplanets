import numpy as np
from scipy.integrate import RK45

from celestial_body import CelestialBody


def propogate_orbits(star: CelestialBody, planets: list[CelestialBody])\
    -> list[CelestialBody]:
    """
    Propogates the orbits of bodies in a planetary system around a single fixed 
    central star.

    Parameters
    ----------
    - star (CelestialBody object): The central body around which the planetary 
    system orbits.
    - planets (list of CelestialBody objects): A list of the planetary bodies 
    that orbit the central star.

    Returns
    -------
    - planets (list of CelestialBody objects): The input planets, updated with 
    new positions and velocities.
    """
