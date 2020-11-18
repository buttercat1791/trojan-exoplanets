from enum import Enum
import numpy as np

from celestial_body import CelestialBody
from propogate_orbits import G


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
        radius: float, position: np.array, velocity: np.array,\
        apoapsis: float, periapsis: float) -> None:
        self.name: str = name
        self.body_type: CelestialType = body_type
        self.mass: float = mass
        self.radius: float = radius
        self.position: np.array = position
        self.velocity: np.array = velocity


    def angular_momentum(self, central: CelestialBody) -> np.array:
        """
        Computes the specific relative angular momentum between the body and 
        the given central body in the system.

        Parameter
        ---------
        - central (CelestialBody): The central body of significantly larger 
        mass around which the body orbits.

        Returns
        -------
        - h (1-d array): The vector of the specific relative angular momentum 
        of the two bodies.
        """

        # Relative position vector.
        r = self.position - central.position
        # Relative velocity vector.
        v = self.velocity - central.velocity
        
        h = np.cross(r, v)

        return h


    def eccentricity(self, central: CelestialBody) -> float:
        """
        Estimates the eccentricity of the body's orbit from it's apoapsis and 
        periapsis relative to the central body.

        Parameter
        ---------
        - central (CelestialBody): The central body of significantly larger 
        mass around which the body orbits.

        Returns
        -------
        - e (float): The eccentricity of the body's orbit around a central body.
        """

        # Specific orbital energy.
        epsilon = self.specific_orbital_energy(central)
        # Standard gravitational parameter.
        mu = self.standard_gravitational_parameter(central)
        # Magnitude of the specific relative angular momentum.
        h = np.linalg.norm(self.angular_momentum(central))

        e = np.sqrt(1 + (2 * epsilon * (h ** 2) / (mu ** 2)))

        return e


    def orbital_period(self, central: CelestialBody) -> float:
        """
        Computes the orbital period of the body relative to a given stationary 
        central body.

        Parameter
        ---------
        - central (CelestialBody): The central body of significantly larger 
        mass around which the body orbits.

        Returns
        -------
        - T (float): The orbital period of the orbiting body.
        """

        # Semi-major axis.
        a = self.semimajor_axis(central)
        # Standard gravitational parameter.
        mu = self.standard_gravitational_parameter(central)

        T = 2 * np.pi * np.sqrt((a ** 3) / mu)

        return mu


    def specific_orbital_energy(self, central: CelestialBody) -> float:
        """
        Computes the specific orbital energy relative to a central body.

        Parameter
        ---------
        - central (CelestialBody): The central body of significantly larger 
        mass around which the body orbits.

        Returns
        -------
        - epsilon (float): The specific orbital energy relative to the passed 
        body.
        """

        # Orbital distance.
        r = np.linalg.norm(self.position - central.position)
        # Relative orbital speed.
        v = np.sqrt(G * central.mass / r)
        # Standard gravitational parameter.
        mu = self.standard_gravitational_parameter(central)

        epsilon = ((v ** 2) / 2) - (mu / r)

        return epsilon


    def semimajor_axis(self, central: CelestialBody) -> float:
        """
        Computes the semi-major axis of the body's orbit around the system's 
        central body.

        Parameter
        ---------
        - central (CelestialBody): The central body of significantly larger 
        mass around which the body orbits.

        Returns
        -------
        - a (float): The approximate semi-major axis relative to the system's 
        central body.
        """

        # Standard gravitational parameter.
        mu = self.standard_gravitational_parameter(central)
        # Specific energy.
        epsilon = self.specific_orbital_energy(central)

        a = -mu / (2 * epsilon)

        return a


    def standard_gravitational_parameter(self, central: CelestialBody) -> float:
        """
        Computes the standard gravitational parameter of the body relative to a 
        central body around which it orbits.

        Parameter
        ---------
        - central (CelestialBody): The central body of significantly larger 
        mass around which the body orbits.

        Returns
        -------
        - mu (float): The standard gravitational parameter of the two masses.
        """

        mu = G * (self.mass + central.mass)

        return mu
