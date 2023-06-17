import MDAnalysis as mda
import numpy as np
import math


def calculate_num_atoms(file, selection):
    """
    Calculates the average number of atoms in a specific region of a trajectory.

    Parameters:
    - file (str): Path to the XYZ file containing the trajectory.
    - selection (str): Selection string for the atoms of interest.

    Returns:
    - num_atoms (float): Average number of atoms in the specified region.

    Requires the MDAnalysis library: https://www.mdanalysis.org/

    The function performs the following steps:
    1. Loads the trajectory and topology.
    2. Selects only the atoms in the specific region.
    3. Calculates the number of atoms in each frame of the trajectory.
    4. Returns the average number of atoms.
    """
    u = mda.Universe(file)

    num_atoms = []
    for ts in u.trajectory[1:]:
        atoms = u.select_atoms(selection)
        num_atoms.append(len(atoms))

    return np.mean(num_atoms)


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


def calculate_metal_area(lattice, atoms):
    """
    Calculates the active metal area from the lattice cell size and the number of metal atoms.

    Parameters:
    - lattice (float): Lattice cell size in Angstrom.
    - atoms (float): Number of metal atoms.

    Returns:
    - metal_area (float): Active metal area in square meters.
    """
    area_atom = lattice ** 2 * math.sqrt(3/4)
    metal_area = (atoms * area_atom) * 1e-20  # m²
    return metal_area


def calculate_total_metal_area(lattice):
    """
    Calculates the total metal area from the lattice cell size.

    Parameters:
    - lattice (float): Lattice cell size in Angstrom.

    Returns:
    - total_area (float): Total metal area in square meters.
    """
    area_atom = lattice ** 2 * math.sqrt(3/4)
    total_area = (492 * area_atom) * 1e-20  # m²
    return total_area


def calculate_permeability(value, area):
    """
    Calculates the permeability from a value and an area.

    Parameters:
    - value (float): Value to consider.
    - area (float): Area in square meters.

    Returns:
    - permeability (float): Permeability in mol/(m²·s).
    """
    permeability = value / area
    return permeability


if __name__ == '__main__':
    file = "/path/to/file.xyz"

    # Calculate the average number of oxygen atoms on the active metal surface
    selection = "name Ox and around 5 name Pt"
    num_atoms = calculate_num_atoms(file, selection)

    # Calculate the number of moles
    n_moles = calculate_n_moles(num_atoms)

    # Calculate the average number of metal atoms
    selection = "name Pt and around 3.0 name O1"
    atoms_Pt = calculate_num_atoms(file, selection)

    # Define the lattice cell size
    lattice_Pt = 3.92  # Angstrom

    # Calculate active and inactive metal areas
    area_active = calculate_metal_area(lattice_Pt, atoms_Pt)
    area_inactive = calculate_metal_area(lattice_Pt, calculate_num_atoms(file, "name Pt and around 3.0 name Cgr C F O3 O4 S"))

    # Calculate total metal area
    total_area = calculate_total_metal_area(lattice_Pt)

    # Calculate metal surface utilization
    utilization = ((total_area - area_inactive) / total_area) * 100

    # Calculate permeability
    value = n_moles / (num_atoms * 10)  # moles per picoseconds
    permeability = calculate_permeability(value, area_active)

    # Generate the results report
    result_report = f"Results Report\n"
    result_report += f"----------------\n"
    result_report += f"Active metal area: {area_active} m²\n"
    result_report += f"Inactive metal area: {area_inactive} m²\n"
    result_report += f"Total metal area: {total_area} m²\n"
    result_report += f"Metal surface utilization: {utilization}%\n"
    result_report += f"Permeability per metal area: {permeability} mol/(m²·s)\n"
    result_report += f"Average number of oxygen atoms on the metal surface in the trajectory: {num_atoms}\n"

    # Save the report to a text file
    with open("result_report.txt", "w") as f:
        f.write(result_report)

    # Print the report to the console
    print(result_report)
