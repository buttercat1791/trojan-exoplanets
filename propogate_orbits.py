import numpy as np
from scipy.integrate import solve_ivp

from celestial_body import CelestialBody


# Newton's Gravitational Constant.
G = 6.67430e-11


def g(index: int, pos: np.array, bodies: list[CelestialBody]):
    """
    Computes the acceleration due to gravity on a given body as a combination 
    of many other celestial bodies.

    Parameters
    ----------
    - index (int): The index of the celestial body whose acceleration is to be 
    computed.
    - pos (1-d array): A 3-d vector of the position of the body at index i.
    - bodies (list of CelestialBody objects): A list of all the bodies in the 
    system, where the first object in the list is the central star.

    Returns
    -------
    - g (1-d array): The acceleration vector of main.
    """

    # Compute the sum total of the accelerations from the other bodies.
    g = np.zeros(3)
    for i in range(1, len(bodies)):
        # Skip main when iterating over the list.
        if i != index:
            # Position vector and mass of the given body.
            mass_body = bodies[i].mass
            pos_body = bodies[i].position
            # Relative position.
            r = pos - pos_body
            # Contribution to main's acceleration.
            g -= G * mass_body * r / (np.linalg.norm(r) ** 3)

    return g


def propogate_orbits(bodies: list[CelestialBody], time_step: float)\
    -> list[CelestialBody]:
    """
    Propogates the orbits of bodies in a planetary system around a single fixed 
    central star.

    Parameters
    ----------
    - bodies (list of CelestialBody objects): A list of all the celestial 
    bodies in the system, where the first object in the list is the central 
    star.
    - time_step (float): The size of the time step in seconds.

    Returns
    -------
    - bodies (list of CelestialBody objects): The input bodies, updated with 
    new positions and velocities.
    """

    # Update the position and velocity of each body in the list.
    for i in range(1, len(bodies)):
        pos = bodies[i].position
        vel = bodies[i].velocity

        # Use RK4 method to compute the new velocity and position.
        k1v = g(i, pos, bodies) * time_step
        k1x = vel * time_step
        k2v = g(i, pos + k1x / 2, bodies) * time_step
        k2x = (vel + k1v / 2) * time_step
        k3v = g(i, pos + k2x / 2) * time_step
        k3x = (vel + k2v / 2) * time_step
        k4v = g(i, pos + k3x, bodies) * time_step
        k4x = (vel + k3v) * time_step
        vel += (k1v + (2 * k2v) + (2 * k3v) + k4v) / 6
        pos += (k1x + (2 * k2x) + (2 * k3x) + k4x) / 6

        # Update the velocity and position.
        bodies[i].position = pos
        bodies[i].velocity = vel

    return bodies
