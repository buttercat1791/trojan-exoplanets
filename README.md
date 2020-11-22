# Trojan Exoplanets
Michael Jurkoic

## Usage
Run the program with the command `Python trojan_exoplanets.py <filename> <step size> <margin>`.
The positional arguments are:
1. `<filename>` The path to a file containing data on the celestial bodies in the system.
2. `<step size>` The size of each simulation time step in seconds. A smaller step increases precision, but makes the code execute slower.
3. `<margin>` A number indicating the percent deviation allowed from a 1:1 resonance between the Trojan pair while running the simulation.

## File Format
The file should be a plaintext file, and may have any number of lines. Each line represents a celestial body, and should take the form:

```
<'STAR' or 'GIANT' or 'TERRESTRIAL'> trojan=<'True' or 'False'> name=<object name> mass=<mass in kilograms> radius=<radius in meters> position=<x>,<y>,<z> velocity=<dx>,<dy>,<dz>
```

Any missing parameters will be set to 0 by default. The `<'STAR' or 'GIANT' or 'TERRESTRIAL'>` parameter must be the first parameter on the line. All other parameters may be listed in any order.
