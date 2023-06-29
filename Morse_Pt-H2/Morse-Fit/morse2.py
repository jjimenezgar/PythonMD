import MDAnalysis as mda
import numpy as np
import pandas as pd

# Cargar datos de simulación de dinámica molecular
u = mda.Universe('/home/johnconnor/Documentos/PtH_MD/morse_parameters/Ea.xyz')

# Seleccionar los átomos para los que se quiere calcular la energía de interacción
atomo1 = u.select_atoms('name H')
atomo2 = u.select_atoms('name Pt1')
atomo3 = u.select_atoms('name Pt2')


n_atomo1=len(atomo1)
n_atomo2=len(atomo2)
n_atomo3=len(atomo3)

 
def morse_2(d_e,alpha,r_e,d_e2,alpha2,r_e2):

    total_energy=[]
    z_distance=[]
    for ts in u.trajectory: ## lee cada frame
            cm = atomo1.center_of_mass() ## Centro de masa
            f_energy=[]
            z_distance.append(cm[2]-11.22)
            for i in range(n_atomo1):
                for x in range(n_atomo2): ### Compara el centro de masa con la superficie
                    r=atomo2[x].position-atomo1[i].position  # Calcula la diferencia vectorial
                    distance=np.linalg.norm(r) #Calcula la distancia   
                    if distance < float(10):
                    # Calcular la energía de interacción utilizando la ecuación de Morse en cada paso de la simulación
                        energias = d_e * (1 - np.exp(-alpha * (distance - r_e)))**2 - d_e

                    # Mostrar los resultados
                        f_energy.append(energias)
                
                for x in range(n_atomo3): ### Compara el centro de masa con la superficie
                    r=atomo3[x].position-atomo1[i].position  # Calcula la diferencia vectorial
                    distance=np.linalg.norm(r) #Calcula la distancia   
                    if distance < float(10):
                    # Calcular la energía de interacción utilizando la ecuación de Morse en cada paso de la simulación
                        energias = d_e2 * (1 - np.exp(-alpha2 * (distance - r_e2)))**2 - d_e2

                    # Mostrar los resultados
                        f_energy.append(energias)
            total_energy.append(sum(f_energy))


    df = pd.DataFrame( { '  ': z_distance,' ': total_energy,})
    
    
    df.to_csv('/home/johnconnor/Documentos/PtH_MD/mp_rep_2/Energies/Energy_{}__{}__{}__{}__{}__{}_.dat'.format(d_e,alpha,r_e,d_e2,alpha2,r_e2),  index=False, sep='\t')


if __name__ == '__main__':

    r_e=2.2

    for d_e in np.arange(0.201, 0.41, 0.01).round(3):
        for d_e2 in [8E-7 ,9E-7 ,1E-6 ,1.1E-6 ,1.2E-6]:   ### 0.8-1.2
            for alpha in np.arange(4.01, 5.2, 0.2).round(3):
                for alpha2 in np.arange(0.301, 1.02, 0.02).round(3):
                    for r_e2 in np.arange(10.01, 20.2, 0.2).round(2):
                        morse_2(d_e,alpha,r_e,d_e2,alpha2,r_e2)


                            
