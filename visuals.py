import numpy as np
from typing import List
import vpython as vp

from celestial_body import CelestialBody, CelestialType


def create_vbodies(bodies: List[CelestialBody]) -> List[vp.sphere]:
    """
    Creates a list of VPython objects from a list of CelestialBody objects.

    Parameters
    ----------
    - bodies (list of CelestialBody objects): The bodies in the system.

    Returns
    -------
    - vbodies (list of VPython objects): A list of VPython objects representing 
    the bodies in the CelestialBody list.
    """

    vbodies = [vp.sphere() for i in len(bodies)]

    # Create the planets.
    for i, body in enumerate(bodies):
        position = vp.vector(body.position[0], body.position[1],\
            body.position[2])
        color = vp.color.blue
        radius = 0
        if body.type == CelestialType.STAR:
            color = vp.color.yellow
            radius = 10e8
        elif body.type == CelestialType.GIANT:
            color = vp.color.red
            radius = 10e7
        elif body.type == CelestialType.TERRESTRIAL:
            color = vp.color.green
            radius = 10e6
        vbodies[i] = vp.sphere(pos=position, radius=radius, color=color,\
            make_trail=True)
    
    return vbodies


def update_vbodies(vbodies: List[vp.sphere], positions: np.array)\
    -> List[vp.sphere]:
    """
    Update the VPython objects to the new positions.

    Parameters
    ----------
    - vbodies (list of VPython objects): The VPython representations of the 
    simulated bodies.
    - positions (2-d array): A list of position vectors corresponding to the 
    bodies in vbodies.

    Returns
    -------
    - vbodies (list of VPython objects): The updated list of VPython 
    representations.
    """

    for i, vbody in enumerate(vbodies):
        new_pos = positions[i]
        pos_vec = vp.vec(new_pos[0], new_pos[1], new_pos[2])
        vbody.pos = pos_vec


def visualize(bodies: List[CelestialBody], positions: np.array):
    """
    Uses the results of a simulation to visualize that simulation in VPython.

    Parameters
    ----------
    - bodies (list of CelestialBody objects): The list of celestial bodies used 
    in the simulation.
    - times (1-d array): The time points that were simulated.
    - positions (3-d array): Contains an array of position vectors 
    corresponding to each time in times.
    """

    # Get a list of VPython objects to represent the simulated bodies.
    vbodies = create_vbodies(bodies)
    # Create a counter for the simulation playback.
    i = 0
    i_max = len(positions)

    # Visualize the simulation.
    while i < i_max:
        vbodies = update_vbodies(vbodies, positions[i])
        i += 1
        vp.rate(100)
