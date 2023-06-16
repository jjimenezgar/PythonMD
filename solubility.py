import numpy as np
import MDAnalysis as mda
import matplotlib.pyplot as plt
import pandas as pd
import free_volume

# Ruta al archivo XYZ
file = "/ruta/al/archivo.xyz"

def calculate_num_atoms(file):
    """
    Calcula el número promedio de átomos de oxígeno en una región específica de una trayectoria.

    Parámetros:
    - file (str): Ruta al archivo XYZ que contiene la trayectoria.

    Retorna:
    - num_oxygen_atoms (int): Número promedio de átomos de oxígeno en la región especificada.

    Requiere la biblioteca MDAnalysis: https://www.mdanalysis.org/

    La función realiza los siguientes pasos:
    1. Carga la trayectoria y la topología.
    2. Selecciona solo los átomos de oxígeno en la región específica.
    3. Calcula el número de átomos de oxígeno en cada paso de la trayectoria.
    4. Retorna el número promedio de átomos de oxígeno.
    """
    # Carga la trayectoria y la topología
    u = mda.Universe(file)

    num_oxygen_atoms = []
    for ts in u.trajectory:
        # Selecciona solo los átomos de oxígeno en la región específica
        oxygen = u.select_atoms("prop z < 30 and name Ox")
        num_oxygen_atoms.append(len(oxygen))

    return np.mean(num_oxygen_atoms)

def calculate_n_moles(n_atoms):
    """
    Calcula el número de moles a partir del número de átomos.

    Parámetros:
    - n_atoms (float): Número de átomos.

    Retorna:
    - n_moles (float): Número de moles.

    Utiliza la constante de Avogadro para realizar la conversión.
    """
    avogadro = 6.022e23
    n_moles = n_atoms / avogadro
    return n_moles

def calculate_concentration(n_moles, volume):
    """
    Calcula la concentración a partir del número de moles y el volumen.

    Parámetros:
    - n_moles (float): Número de moles.
    - volume (float): Volumen en metros cúbicos.

    Retorna:
    - concentration (float): Concentración en moles por metro cúbico.
    """
    concentration = n_moles / volume
    return concentration

def calculate_pressure(n_moles, volume):
    """
    Calcula la presión a partir del número de moles y el volumen.

    Parámetros:
    - n_moles (float): Número de moles.
    - volume (float): Volumen en metros cúbicos.

    Retorna:
    - pressure (float): Presión en atmósferas.
    """
    R = 8.314
    temperature = 298.15  # Kelvin
    pressure_pa = (3.98e-21 * R * temperature) / 1.35e-24  # nmol para 2400 partículas y el volumen es el de la caja menos la esfera
    pressure_atm = pressure_pa / 101325.0  # Convertir Pascal a atmósferas
    return pressure_atm

def calculate_solubility(concentration, pressure):
    """
    Calcula la solubilidad a partir de la concentración y la presión.

    Parámetros:
    - concentration (float): Concentración en moles por metro cúbico.
    - pressure (float): Presión en atmósferas.

    Retorna:
    - solubility (float): Solubilidad en moles por metro cúbico por atmósfera.
    """
    solubility = concentration / pressure
    return solubility

# Calcular propiedades
volume = free_volume.calculate_free_volume(file, 30)
n_atoms = calculate_num_atoms(file)
n_moles = calculate_n_moles(n_atoms)
C = calculate_concentration(n_moles, volume[1])  # moles por metro cúbico
P = calculate_pressure(n_moles, volume[1])  # atm
S = calculate_solubility(C, P)  # moles por metro cúbico por atm

# Imprimir el reporte
print("Reporte de Propiedades")
print("----------------------")
print("Volumen libre: ", volume[0], " m³")
print("Volumen total: ", volume[1], " m³")
print("Número de átomos: ", n_atoms)
print("Número de moles: ", n_moles, " mol")
print("Concentración: ", C, " mol/m³")
print("Presión: ", P, " atm")
print("Solubilidad: ", S, " mol/(m³·atm)")
