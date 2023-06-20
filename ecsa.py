import numpy as np
import MDAnalysis as mda

class ECSAnalysis:
    def __init__(self, file1, file2,Lattice,molar_mass):
        self.file1 = file1
        self.file2 = file2
        self.Lattice = float(Lattice)  # Angstrom
        self.molar_mass = float(molar_mass)  # g/mol of Pt
        self.Atomic_area_CN9 = (self.Lattice ** 2 * np.sqrt(3)) / 4  # Atomic area for Pt (111) surface
        self.Atomic_area_CN8 = (self.Lattice ** 2) / 2  # Atomic area for Pt (100) surface
        self.Atomic_area_CN6 = (self.Lattice ** 2 * np.sqrt(3)) / 2  # Atomic area for Pt NP vertices

    def surface_contact(self, atom1, atom2, r_cut):
        """
        Calculates the number of contacts between two sets of atoms in different files.

        Args:
            atom1 (str): Atom selection string for the first set of atoms.
            atom2 (str): Atom selection string for the second set of atoms.
            r_cut (float): Cutoff distance for considering two atoms in contact.

        Returns:
            float: Average number of contacts between the two sets of atoms.
        """
        u1 = mda.Universe(self.file1)
        u2 = mda.Universe(self.file2)

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

    def SA_metal(self):
        """
        Calculates the surface area of a nanoparticle and the total number of atoms in the NP.

        Returns:
            tuple: A tuple containing the total surface area of the nanoparticle and the total number of atoms in the NP.

        Note:
            The function assumes a specific lattice parameter and atomic areas for different surface types.

        References:
            - Lattice parameter and atomic areas: [provide reference/source]
        """
        u1 = mda.Universe(self.file1)

        total_atoms = u1.select_atoms('name 6 8 9 12')  # Select surface atoms with coordination number 6, 8, 9, 12
        CN_6 = u1.select_atoms('name 6')  # Select surface atoms with coordination number 6
        CN_8 = u1.select_atoms('name 8')  # Select surface atoms with coordination number 8
        CN_9 = u1.select_atoms('name 9')  # Select surface atoms with coordination number 9

        number_atom6 = len(CN_6)
        number_atom8 = len(CN_8)
        number_atom9 = len(CN_9)
        number_total = len(total_atoms)

        SA_metal = (
            number_atom6 * self.Atomic_area_CN6 +
            number_atom8 * self.Atomic_area_CN8 +
            number_atom9 * self.Atomic_area_CN9
        ) * 1e-20  # Total nanoparticle area

        return SA_metal, number_total

    def calculate_surface_area(self, atom2, r_cut):
        """
        Calculates the surface area of a nanoparticle and the number of contacts for different types of atoms.

        Args:
            atom2 (str): Atom selection string for the second set of atoms.
            r_cut (float): Cutoff distance for considering two atoms in contact.

        Returns:
            tuple: A tuple containing the total surface area, the surface area of the nanoparticle,
                   and the total number of contacts for different types of atoms.
        """
        CN = ["6", "8", "9"]
        area_total = 0
        number_contact = 0

        for i in CN:
            numero = self.surface_contact(i, atom2, r_cut)

            if i == "6":
                SA = numero * self.Atomic_area_CN6 * 1e-20  # Convert Å^2 to m^2
            elif i == "8":
                SA = numero * self.Atomic_area_CN8 * 1e-20  # Convert Å^2 to m^2
            elif i == "9":
                SA = numero * self.Atomic_area_CN9 * 1e-20  # Convert Å^2 to m^2

            area_total += SA
            number_contact += numero

        return area_total, number_contact

    def ecsa_maxima(self, surface_NP, total_atoms):
        """
        Calculates the maximum electrochemical surface area (ECSA) for a nanoparticle.

        Args:
            surface_NP (float): Surface area of the nanoparticle.

        Returns:
            float: Maximum ECSA in m^2/g.
        """
        ma_metal = self.molar_mass  # g/mol of Pt
        NAvo = 6.022e23  # mol^-1
        n_atoms = total_atoms  # Number of atoms in the icosahedral nanoparticle
        ecsa_maxima = surface_NP / (n_atoms * (ma_metal / NAvo))  # m^2/g
        return ecsa_maxima

    def ecsa_active(self, active_surface, ecsa_maxima, total_atoms):
        """
        Calculates the active electrochemical surface area (ECSA) and utilization percentage.

        Args:
            active_surface (float): Active surface area.
            ecsa_maxima (float): Maximum ECSA.

        Returns:
            tuple: A tuple containing the active ECSA and the utilization percentage.
        """
        ma_metal = self.molar_mass  # g/mol of Pt
        NAvo = 6.022e23  # mol^-1
        n_atoms = total_atoms  # Number of atoms in the icosahedral nanoparticle
        ecsa_active = active_surface / (n_atoms * (ma_metal / NAvo))  # m^2/g
        utilization = (ecsa_active / ecsa_maxima) * 100  # Catalyst utilization percentage
        return ecsa_active, utilization

if __name__ == '__main__':
    file1 = input("Enter the path to the first file (Metal_CN): ")
    file2 = input("Enter the path to the second file(Rest of the system): ")
    Lattice = input("Enter the lattice parameter of the nanoparticle: ")  ## 3.92 Å for Pt
    molar_mass = input("Enter the molar mass of the metal: ")  ## 195 g/mol for Pt
    atom2 = input("Enter the atom related to the active area: ")
    r_cut = float(input("Enter the cutoff distance for considering two atoms in contact: "))


    ECSA = ECSAnalysis(file1, file2,Lattice,molar_mass) ## Create an instance of the NanoparticleAnalysis class.
    Surface_metal, total_atoms = ECSA.SA_metal() ## Get the total surface area of the nanoparticle and the total number of atoms in the NP.
    ecsa_maxima = ECSA.ecsa_maxima(Surface_metal, total_atoms) ## Get the ecsa_maxima
    active_surface,n_atoms = ECSA.calculate_surface_area(atom2,r_cut) ## Get the active surface area and the total number of atoms in the active area.
    ecsa_active, utilization = ECSA.ecsa_active(active_surface, ecsa_maxima, total_atoms) ## Get the ecsa_active and the utilization percentage.

    print("ECSA Active: ", ecsa_active, "m²/g")
    print("Utilization: ", utilization, "%")
