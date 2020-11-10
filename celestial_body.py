from enum import Enum
import numpy as np


class CelestialType(Enum):
    STAR = 0
    GIANT = 1
    TERRESTRIAL = 2


class CelestialBody:
    """
    Class representation of a celestial body.

    Constructor Parameters
    ----------------------
    - name (string): The name of the celestial body being represented.
    - body_type (CelestialType): An enum member representing the type of 
    celestial body being represented.
    - mass (float): The mass in kilograms of the body being represented.
    - radius (float): The radius in meters of the body being represented.
    - position (1-d array): A 3-d vector of the form [x, y, z] representing 
    the initial position of the body in meters from the origin along the x, 
    y, and z axes.
    - velocity (1-d array): A 3-d vector of the form [vx, vy, vz] 
    representing the initial velocity of the body in the x, y, and z 
    directions.
    """

    def __init__(self, name: str, body_type: CelestialType, mass: float,\
        radius: float, position: np.array, velocity: np.array) -> None:
        self.name: str = name
        self.body_type: CelestialType = body_type
        self.mass: float = mass
        self.radius: float = radius
        self.position: np.array = position
        self.velocity: np.array = velocity
