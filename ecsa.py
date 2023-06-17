import numpy as np
import MDAnalysis as mda


def surface_contact(file1, file2, atom1, atom2, r_cut):
    """
    Calculates the number of contacts between two sets of atoms in different files.

    Args:
        file1 (str): Path to the first file.
        file2 (str): Path to the second file.
        atom1 (str): Atom selection string for the first set of atoms.
        atom2 (str): Atom selection string for the second set of atoms.
        r_cut (float): Cutoff distance for considering two atoms in contact.

    Returns:
        float: Average number of contacts between the two sets of atoms.
    """
    u1 = mda.Universe(file1)
    u2 = mda.Universe(file2)

    term1 = u1.select_atoms(f'name {atom1}')  # Select surface atoms
    term2 = u2.select_atoms(f'name {atom2}')  # Select atoms to calculate center of mass

    number_atom1 = len(term1)
    number_atom2 = len(term2)

    list_contact = []
    for ts in u1.trajectory:
        if ts.frame == 0:
            count = 0
            for x in range(number_atom1):
                list_d = []
                for i in range(number_atom2):
                    r = term1[x].position - term2[i].position
                    d = np.linalg.norm(r)
                    if d < r_cut:
                        list_d.append(d)

                if len(list_d) > 0 and min(list_d) < r_cut:
                    count += 1

            list_contact.append(count)

    return np.mean(list_contact)


def calculate_surface_area(file1, file2, atom2, r_cut):
    """
    Calculates the surface area of a nanoparticle and the number of contacts for different types of atoms.

    Args:
        file1 (str): Path to the first file.
        file2 (str): Path to the second file.
        atom2 (str): Atom selection string for the second set of atoms.
        r_cut (float): Cutoff distance for considering two atoms in contact.

    Returns:
        tuple: A tuple containing the total surface area, the surface area of the nanoparticle,
               and the total number of contacts for different types of atoms.
    """
    Lattice_Pt = 3.92  # Angstrom
    Atomic_area_CN9 = (Lattice_Pt ** 2 * np.sqrt(3)) / 4  # Atomic area for Pt (111) surface
    Atomic_area_CN8 = (Lattice_Pt ** 2) / 2  # Atomic area for Pt (100) surface
    Atomic_area_CN6 = (Lattice_Pt ** 2 * np.sqrt(3)) / 2  # Atomic area for Pt NP vertices

    SA_1415 = (12 * Atomic_area_CN6 + 184 * Atomic_area_CN8 + 296 * Atomic_area_CN9) * 1e-20  # Total nanoparticle area

    CN = ["6", "8", "9"]
    area_total = 0
    number_contact = 0
    for i in CN:
        numero = surface_contact(file1, file2, i, atom2, r_cut)
        if i == "6":
            SA = numero * Atomic_area_CN6 * 1e-20  # Convert Å^2 to m^2
        elif i == "8":
            SA = numero * Atomic_area_CN8 * 1e-20  # Convert Å^2 to m^2
        elif i == "9":
            SA = numero * Atomic_area_CN9 * 1e-20  # Convert Å^2 to m^2
        area_total += SA
        number_contact += numero

    return area_total, SA_1415, number_contact


def ecsa_maxima(surface_NP):
    """
    Calculates the maximum electrochemical surface area (ECSA) for a nanoparticle.

    Args:
        surface_NP (float): Surface area of the nanoparticle.

    Returns:
        float: Maximum ECSA in m^2/g.
    """
    ma_metal = 195  # g/mol of Pt
    NAvo = 6.022e23  # mol^-1
    n_atoms = 1415  # Number of atoms in the icosahedral nanoparticle
    ecsa_maxima = surface_NP / (n_atoms * (ma_metal / NAvo))  # m^2/g
    return ecsa_maxima


def ecsa_active(active_surface, ecsa_maxima):
    """
    Calculates the active electrochemical surface area (ECSA) and utilization percentage.

    Args:
        active_surface (float): Active surface area.
        ecsa_maxima (float): Maximum ECSA.

    Returns:
        tuple: A tuple containing the active ECSA and the utilization percentage.
    """
    ma_metal = 195  # g/mol of Pt
    NAvo = 6.022e23  # mol^-1
    n_atoms = 1415  # Number of atoms in the icosahedral nanoparticle
    ecsa_active = active_surface / (n_atoms * (ma_metal / NAvo))  # m^2/g
    utilization = (ecsa_active / ecsa_maxima) * 100  # Catalyst utilization percentage
    return ecsa_active, utilization


if __name__ == '__main__':
    file1 = "/home/ruta/archivo1.xyz"  ## archivo1.xyz is the file with the surface metal atoms, identify for coordination number
    file2 = "/home/ruta/archivo2.xyz" ## archivo2.xyz is the file with the rest of the atoms of the system (water, carbon vulcan,Nafion, etc)
    
    atom2 = "O1"  ## Atom related to active area (water molecule)
    r_cut = 3.2 ## Cutoff distance for considering two atoms in contact

    
    total_area = calculate_surface_area(file1, file2, atom2, r_cut)  ## Area total, area nanoparticle, number of contacts

    surface_maxima = ecsa_maxima(total_area[1]) ## ECSA maxima

    ecsa_active_result = ecsa_active(total_area[0], surface_maxima) ## ECSA active and utilization

    ecsa_active_surface = ecsa_active_result[0]  ## ECSA active
    ecsa_active_utilization = ecsa_active_result[1] # ECSA utilization

    print(f"ECSA Active Surface: {ecsa_active_surface:.4f} m^2/g")
    print(f"ECSA Active Utilization: {ecsa_active_utilization:.2f}%")
