import MDAnalysis as mda
import numpy as np
import math


def calculate_num_atoms(file, selection):
    """
    Calcula el número promedio de átomos en una región específica de una trayectoria.

    Parámetros:
    - file (str): Ruta al archivo XYZ que contiene la trayectoria.
    - selection (str): Cadena de selección para los átomos de interés.

    Retorna:
    - num_atoms (float): Número promedio de átomos en la región especificada.

    Requiere la biblioteca MDAnalysis: https://www.mdanalysis.org/

    La función realiza los siguientes pasos:
    1. Carga la trayectoria y la topología.
    2. Selecciona solo los átomos en la región específica.
    3. Calcula el número de átomos en cada paso de la trayectoria.
    4. Retorna el número promedio de átomos.
    """
    u = mda.Universe(file)

    num_atoms = []
    for ts in u.trajectory[1:]:
        atoms = u.select_atoms(selection)
        num_atoms.append(len(atoms))

    return np.mean(num_atoms)


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


def calculate_metal_area(lattice, atoms):
    """
    Calcula el área metálica activa a partir del tamaño de la celda de la red y el número de átomos metálicos.

    Parámetros:
    - lattice (float): Tamaño de la celda de la red en Angstrom.
    - atoms (float): Número de átomos metálicos.

    Retorna:
    - metal_area (float): Área metálica activa en metros cuadrados.
    """
    area_atom = lattice ** 2 * math.sqrt(3/4)
    metal_area = (atoms * area_atom) * 1e-20  # m²
    return metal_area


def calculate_total_metal_area(lattice):
    """
    Calcula el área metálica total a partir del tamaño de la celda de la red.

    Parámetros:
    - lattice (float): Tamaño de la celda de la red en Angstrom.

    Retorna:
    - total_area (float): Área metálica total en metros cuadrados.
    """
    area_atom = lattice ** 2 * math.sqrt(3/4)
    total_area = (492 * area_atom) * 1e-20  # m²
    return total_area


def calculate_permeability(value, area):
    """
    Calcula la permeabilidad a partir de un valor y un área.

    Parámetros:
    - value (float): Valor a considerar.
    - area (float): Área en metros cuadrados.

    Retorna:
    - permeability (float): Permeabilidad en mol/(m²·s).
    """
    permeability = value / area
    return permeability


file = "/ruta/al/archivo.xyz"

# Calcular número promedio de átomos de oxígeno en la superficie metálica activa
selection = "name Ox and around 5 name Pt"
num_atoms = calculate_num_atoms(file, selection)

# Calcular número de moles
n_moles = calculate_n_moles(num_atoms)

# Calcular número promedio de átomos metálicos
selection = "name Pt and around 3.0 name O1"
atoms_Pt = calculate_num_atoms(file, selection)

# Definir tamaño de la celda de la red
lattice_Pt = 3.92  # Angstrom

# Calcular área metálica activa e inactiva
area_activa = calculate_metal_area(lattice_Pt, atoms_Pt)
area_inactiva = calculate_metal_area(lattice_Pt, calculate_num_atoms(file, "name Pt and around 3.0 name Cgr C F O3 O4 S"))

# Calcular área metálica total
total_area = calculate_total_metal_area(lattice_Pt)

# Calcular utilización de la superficie metálica
utilization = ((total_area - area_inactiva) / total_area) * 100

# Calcular permeabilidad
value = n_moles / (num_atoms * 10)  # moles por picosegundos
permeability = calculate_permeability(value, area_activa)

# Generar el reporte de resultados
reporte_resultados = f"Reporte de Resultados\n"
reporte_resultados += f"---------------------\n"
reporte_resultados += f"Area metálica activa: {area_activa} m²\n"
reporte_resultados += f"Area metálica inactiva: {area_inactiva} m²\n"
reporte_resultados += f"Area metálica total: {total_area} m²\n"
reporte_resultados += f"Utilización de la superficie metálica: {utilization}%\n"
reporte_resultados += f"Permeabilidad por área metálica: {permeability} mol/(m²·s)\n"
reporte_resultados += f"Número promedio de átomos de oxígeno en la superficie metálica en la dinámica: {num_atoms}\n"

# Guardar el reporte en un archivo de texto
with open("reporte_resultados.txt", "w") as f:
    f.write(reporte_resultados)

# Imprimir el reporte en la consola
print(reporte_resultados)
