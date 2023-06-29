import MDAnalysis as mda
import numpy as np
import pandas as pd

def morse_2(d_e, alpha, r_e, d_e2, alpha2, r_e2):
    """
    Calculate interaction energies using the Morse potential for specified atom selections.

    Parameters:
    - d_e (float): Dissociation energy parameter for the first interaction.
    - alpha (float): Force constant parameter for the first interaction.
    - r_e (float): Equilibrium distance parameter for the first interaction.
    - d_e2 (float): Dissociation energy parameter for the second interaction.
    - alpha2 (float): Force constant parameter for the second interaction.
    - r_e2 (float): Equilibrium distance parameter for the second interaction.
    """

    # Load molecular dynamics simulation data
    u = mda.Universe('/home/johnconnor/Documentos/PtH_MD/morse_parameters/Ea.xyz')

    # Select atoms for which interaction energy will be calculated
    atomo1 = u.select_atoms('name H')
    atomo2 = u.select_atoms('name Pt1')
    atomo3 = u.select_atoms('name Pt2')

    n_atomo1 = len(atomo1)
    n_atomo2 = len(atomo2)
    n_atomo3 = len(atomo3)

    total_energy = []
    z_distance = []

    for ts in u.trajectory:
        cm = atomo1.center_of_mass()

        f_energy = []
        z_distance.append(cm[2] - 11.22)

        for i in range(n_atomo1):
            for x in range(n_atomo2):
                r = atomo2[x].position - atomo1[i].position
                distance = np.linalg.norm(r)
                
                if distance < float(10):
                    energias = d_e * (1 - np.exp(-alpha * (distance - r_e)))**2 - d_e
                    f_energy.append(energias)
                
            for x in range(n_atomo3):
                r = atomo3[x].position - atomo1[i].position
                distance = np.linalg.norm(r)
                
                if distance < float(10):
                    energias = d_e2 * (1 - np.exp(-alpha2 * (distance - r_e2)))**2 - d_e2
                    f_energy.append(energias)
        
        total_energy.append(sum(f_energy))

    df = pd.DataFrame({'Distance': z_distance, 'Energy': total_energy})
    df.to_csv('/home/johnconnor/Documentos/PtH_MD/mp_rep_2/Energies/Energy_{}__{}__{}__{}__{}__{}_.dat'.format(d_e, alpha, r_e, d_e2, alpha2, r_e2), index=False, sep='\t')

if __name__ == '__main__':
    r_e = 2.2

    for d_e in np.arange(0.201, 0.41, 0.01).round(3):
        for d_e2 in [8E-7, 9E-7, 1E-6, 1.1E-6, 1.2E-6]:
            for alpha in np.arange(4.01, 5.2, 0.2).round(3):
                for alpha2 in np.arange(0.301, 1.02, 0.02).round(3):
                    for r_e2 in np.arange(10.01, 20.2, 0.2).round(2):
                        morse_2(d_e, alpha, r_e, d_e2, alpha2, r_e2)


                            
