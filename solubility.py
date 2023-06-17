import numpy as np
import MDAnalysis as mda
import matplotlib.pyplot as plt
import pandas as pd
import free_volume

# Path to the XYZ file
file = "/path/to/file.xyz"

def calculate_num_atoms(file):
    """
    Calculates the average number of oxygen atoms in a specific region of a trajectory.

    Parameters:
    - file (str): Path to the XYZ file containing the trajectory.

    Returns:
    - num_oxygen_atoms (int): Average number of oxygen atoms in the specified region.

    Requires the MDAnalysis library: https://www.mdanalysis.org/

    The function performs the following steps:
    1. Load the trajectory and topology.
    2. Select only the oxygen atoms in the specific region.
    3. Calculate the number of oxygen atoms in each frame of the trajectory.
    4. Return the average number of oxygen atoms.
    """
    # Load the trajectory and topology
    u = mda.Universe(file)

    num_oxygen_atoms = []
    for ts in u.trajectory:
        # Select only the oxygen atoms in the specific region
        oxygen = u.select_atoms("prop z < 30 and name Ox")
        num_oxygen_atoms.append(len(oxygen))

    return np.mean(num_oxygen_atoms)

def calculate_n_moles(n_atoms):
    """
    Calculates the number of moles from the number of atoms.

    Parameters:
    - n_atoms (float): Number of atoms.

    Returns:
    - n_moles (float): Number of moles.

    Uses Avogadro's constant for the conversion.
    """
    avogadro = 6.022e23
    n_moles = n_atoms / avogadro
    return n_moles

def calculate_concentration(n_moles, volume):
    """
    Calculates the concentration from the number of moles and volume.

    Parameters:
    - n_moles (float): Number of moles.
    - volume (float): Volume in cubic meters.

    Returns:
    - concentration (float): Concentration in moles per cubic meter.
    """
    concentration = n_moles / volume
    return concentration

def calculate_pressure(n_moles, volume):
    """
    Calculates the pressure from the number of moles and volume.

    Parameters:
    - n_moles (float): Number of moles.
    - volume (float): Volume in cubic meters.

    Returns:
    - pressure (float): Pressure in atmospheres.
    """
    R = 8.314
    temperature = 298.15  # Kelvin
    pressure_pa = (3.98e-21 * R * temperature) / 1.35e-24  # nmol for 2400 particles and the volume is the box minus the sphere
    pressure_atm = pressure_pa / 101325.0  # Convert Pascal to atmospheres
    return pressure_atm

def calculate_solubility(concentration, pressure):
    """
    Calculates the solubility from the concentration and pressure.

    Parameters:
    - concentration (float): Concentration in moles per cubic meter.
    - pressure (float): Pressure in atmospheres.

    Returns:
    - solubility (float): Solubility in moles per cubic meter per atmosphere.
    """
    solubility = concentration / pressure
    return solubility


if __name__ == "__main__":
    # Calculate properties
    volume = free_volume.calculate_free_volume(file, 30)
    n_atoms = calculate_num_atoms(file)
    n_moles = calculate_n_moles(n_atoms)
    C = calculate_concentration(n_moles, volume[1])  # moles per cubic meter
    P = calculate_pressure(n_moles, volume[1])  # atm
    S = calculate_solubility(C, P)  # moles per cubic meter per atm

    # Print the report
    print("Properties Report")
    print("-----------------")
    print("Free Volume: ", volume[0], " m³")
    print("Total Volume: ", volume[1], " m³")
    print("Number of Atoms: ", n_atoms)
    print("Number of Moles: ", n_moles, " mol")
    print("Concentration: ", C, " mol/m³")
    print("Pressure: ", P, " atm")
    print("Solubility: ", S, " mol/(m³·atm)")
