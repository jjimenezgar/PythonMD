import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis import density
import matplotlib.pyplot as plt

# Archivo de entrada
file = "input.xyz"

# Crear objeto de Universe
u = mda.Universe(file)

def density_2D(u):
    """
    Calcula y visualiza la densidad media 2D de agua en una caja de simulación.

    Args:
        u (Universe): Objeto Universe de MDAnalysis.

    Returns:
        None
    """

    # Definir la coordenada Z máxima de la caja
    z_hi = "Introducir maximo de caja en Z"

    # Seleccionar los átomos de oxígeno en la región deseada (por encima de z_hi)
    ow = u.select_atoms("name O1 and prop {} > z".format(z_hi))

    # Realizar el análisis de densidad
    dens = density.DensityAnalysis(ow, delta=3.5, padding=0)
    dens.run()

    # Obtener la grilla de densidad y convertirla a unidades SPC/E
    grid = dens.results.density.grid
    dens.results.density.convert_density('SPC')

    # Calcular la densidad media en el eje Z (vista superior)
    avg = grid.mean(axis=-1)

    # Crear la figura y los ejes
    fig, ax = plt.subplots(figsize=(10, 8))

    # Crear el mapa de calor de la densidad media
    im = ax.imshow(avg, interpolation="bicubic", cmap='jet', vmin=0, vmax=1)

    # Agregar una barra de color
    cbar = plt.colorbar(im)
    cbar.set_label('Mean density of water over SPC/E literature value')

    # Etiquetas de los ejes
    plt.xlabel('X-axis ($\AA$)')
    plt.ylabel('Y-axis ($\AA$)')

# Llamar a la función para ejecutar el análisis
density_2D(u)

