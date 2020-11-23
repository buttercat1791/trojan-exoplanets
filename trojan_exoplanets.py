import argparse
import matplotlib.pyplot as plt
import numpy as np
from typing import List, Union

from celestial_body import CelestialBody, CelestialType
from propogate_orbits import g, propogate_Verlet


def parse_line(line: str) -> CelestialBody:
    """
    Parses a single line of text to produce a CelestialBody object.

    Parameter
    ---------
    - line (string): A line of text containing the information necessary to 
    construct a CelestialBody object.

    Returns
    -------
    - body (CelestialBody): The CelestialBody object constructed from the 
    parameters given by the input line.
    """

    params: List[str] = line.split(" ")
    
    # Variables used to construct the CelestialBody object.
    body_type: CelestialType = -1
    body_is_trojan: bool = False
    body_name: str = ""
    body_mass: float = 0.0
    body_radius: float = 0.0
    body_position: np.array = np.zeros(3)
    body_velocity: np.array = np.zeros(3)

    # The parsing function cannot detect missing parameters.
    for i, param in enumerate(params):
        if i == 0:
            # Each line should start by specifying the body type.
            if params[0] == "STAR":
                body_type = CelestialType.STAR
            elif params[0] == "GIANT":
                body_type = CelestialType.GIANT
            elif params[0] == "TERRESTRIAL":
                body_type = CelestialType.TERRESTRIAL
            else:
                raise Exception("Each line must start with a valid object type"\
                    + " specification.")
        else:
            tokens = param.split("=")
            if tokens[0] == "trojan":
                body_is_trojan = tokens[1] == "True"
            elif tokens[0] == "name":
                body_name = tokens[1]
            elif tokens[0] == "mass":
                body_mass = float(tokens[1])
            elif tokens[0] == "radius":
                body_radius = float(tokens[1])
            elif tokens[0] == "position":
                vals = tokens[1].split(",")
                if len(vals) != 3:
                    raise Exception("Position must be specified by three"\
                        + " comma-separated floats")
                body_position = np.array([float(vals[0]), float(vals[1]),\
                    float(vals[2])])
            elif tokens[0] == "velocity":
                vals = tokens[1].split(",")
                if len(vals) != 3:
                    raise Exception("Velocity must be specified by three"\
                        + " comma-separated floats")
                body_velocity = np.array([float(vals[0]), float(vals[1]),\
                    float(vals[2])])
            else:
                raise Exception("Invalid token detected while parsing input"\
                    + " file.")

    body = CelestialBody(type=body_type, trojan=body_is_trojan, name=body_name,\
        mass=body_mass, radius=body_radius, position=body_position,\
        velocity=body_velocity)

    return body
    

def parse_system(file: str) -> Union[list, list]:
    """
    Initializes the starting conditions of the simulation from a file with the 
    given name.

    Parameter
    ---------
    - file (string): The name of the file to parse from.

    Returns
    -------
    - bodies (list of CelestialBody objects): A list of celestial bodies in the 
    system, initialized with position, velocity, mass, and radius.
    - trojans (2-tuple of ints): The list indices of the Trojan pair.
    """

    infile = open(file, "r")

    bodies: List[CelestialBody] = list()

    lines: List[str] = infile.readlines()
    star_count = 0
    trojan_count = 0
    trojans = [0, 0]
    for i, line in enumerate(lines):
        body = parse_line(line)
        # Keep track of the number of stars in the system.
        if body.type == CelestialType.STAR:
            star_count += 1
        # Keep track of the number of trojans in the system, and store the 
        # indices of the Trojan pair.
        if body.trojan:
            trojans[trojan_count] = i
            trojan_count += 1
        # We only want to work with single-star systems.
        if star_count > 1:
            raise Exception("Multi-star systems are not allowed.")
        bodies.append(body)
    # We only want to have one Trojan pair at a time.
    if trojan_count != 2:
        raise Exception("The system must have a single Trojan pair.")

    # Put the star in the front of the list, if necessary.
    if bodies[0].type != CelestialType.STAR:
        for i in range(len(bodies)):
            if bodies[i].type == CelestialType.STAR:
                bodies[0], bodies[i] = bodies[i], bodies[0]

    return bodies, trojans


def check_margins(trojan1: CelestialBody, trojan2: CelestialBody,\
    central: CelestialBody, margin: float) -> bool:
    """
    Checks to see if two planets in a Trojan pair are within a given percent 
    margin of a 1:1 period resonance.

    Parameters
    ----------
    - trojan1 (CelestialBody): One member of the co-orbital pair.
    - trojan2 (CelestialBody): The other member of the co-orbital pair.
    - central (CelestialBody): The central body in the system.
    - margin (float): The allowed percent deviation from a 1:1 period resonance.

    Returns
    -------
    - in_margin (bool): True if within the margin, false if not.
    """

    # Get the orbital periods of each planet.
    P1 = trojan1.period(central)
    P2 = trojan2.period(central)

    # Compute the percent difference between the periods.
    diff = abs(P1 - P2)
    avg = np.average([P1, P2])
    p_diff = (diff / avg) * 100

    # Will return true if the percent difference exceeds the margin.
    return not (p_diff > margin)


def plot_periods(t: np.array, p1: np.array, p1_name: str, p2: np.array,\
    p2_name: str):
    """
    Plots the periods vs. time of the given Trojan pair.

    Parameters
    ----------
    - t (1-d array): An array of time values, where each value represents a 
    number of years that have passed.
    - p1 (1-d array): An array of the same length as t giving the orbital 
    period in days of one member of the Trojan pair for each year in t.
    - p1_name (string): The name of one member of the Trojan pair.
    - p2 (1-d array): An array of the same length as t giving the orbital 
    period in days of the other member of the Trojan pair for each year in t.
    - p2_name (string): The name of the other member of the Trojan pair.
    """

    plt.title("Periods vs. Time")
    plt.xlabel("Time (Years)")
    plt.ylabel("Period (Days)")
    plt.plot(t, p1, 'r-', label=p1_name)
    plt.plot(t, p2, 'b-', label=p2_name)
    plt.legend()
    plt.show()


def simulate(bodies: List[CelestialBody], time_step: int,\
    trojans: list, margin: float):
    """
    Takes the parsed parameters and runs a simulation with them.

    Parameters
    ----------
    - bodies (list of CelestialBody objects): A list containing all the bodies 
    in the system with their initial positions and velocities. The system's 
    central star should be at position 0.
    - time_step (int): The size of the simulation time step in seconds. The 
    shorter the time step, the more precise the simulation will be, but the 
    longer it will take to complete.
    - trojans (2-tuple of ints): The indices within the bodies list of the 
    co-orbital exoplanets being investigated.
    - margin (float): The percent deviation allowed from a 1:1 resonance 
    between the Trojan planets before the simulation is stopped.
    """

    # Length of a year in seconds.
    YEAR: float = 60 * 60 * 24 * 365.25
    # Length of a day in seconds.
    DAY = 60 * 60 * 24
    # Time counter (holds time elapsed in seconds).
    time: float = 0
    # Year counter (counts number of years elapsed).
    years: int = 0
    # Set to false to end the simulation.
    in_margin: bool = True

    # Arrays to hold the elapsed number of years and the corresponding orbital 
    # periods to track change over time.
    t = np.array([])
    p1 = np.array([])
    p2 = np.array([])

    # Initialize the accelerations.
    for i, body in enumerate(bodies):
        body.acceleration = g(i, body.position, bodies)

    while in_margin:
        # Move the planets (the star stays still) and update the time.
        bodies = propogate_Verlet(bodies=bodies, time_step=time_step)
        time += time_step

        # Round down to get the number of elapsed years. Update the years 
        # counter if necessary.
        if int(time / YEAR) > years:
            years += 1
            # Update our output arrays appropriately.
            t = np.append(t, years)
            p1 = np.append(p1, bodies[trojans[0]].period(bodies[0]) / DAY)
            p2 = np.append(p2, bodies[trojans[1]].period(bodies[0]) / DAY)
            # Check once a year to see if the Trojan pair is within margins.
            in_margin = check_margins(bodies[trojans[0]], bodies[trojans[1]],\
                bodies[0], margin)
            # Print a status to indicate the program is working.
            if years % 1000 == 0:
                print(f"{years} years elapsed")
                print(f"P1 = {bodies[trojans[0]].period(bodies[0]) / DAY}")
                print(f"P2 = {bodies[trojans[1]].period(bodies[0]) / DAY}")

        # End the loop after a million years.
        if years >= 10e6:
            in_margin = False

    # After the simulation loop finishes, indicate how long it lasted.
    print(f"\nThe Trojan pair remained stable for {years} years.")

    # Plot the change in orbital periods over time.
    plot_periods(t, p1, bodies[trojans[0]].name, p2, bodies[trojans[1]].name)


if __name__ == "__main__":
    # Parse command line arguments.
    parser = argparse.ArgumentParser()
    parser.add_argument("file",\
        help="Extract starting conditions from the given file", type=str)
    parser.add_argument("step",\
        help="Specify the simulation time step size in seconds", type=int)
    parser.add_argument("margin", help="Specify the allowed percent deviation"\
        + " from a 1:1 resonance in the Trojan pair", type=float)
    args = parser.parse_args()

    # Construct the system from the file, and run the simulation.
    bodies, trojans = parse_system(args.file)
    simulate(bodies, args.step, trojans, args.margin)
