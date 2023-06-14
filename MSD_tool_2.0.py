
import MDAnalysis.analysis.msd as msda
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
import seaborn as sns
import pandas as pd

def msd_z(file, atom, z_lo, z_hi):
    """
    Calcula el desplazamiento cuadrático medio (MSD) de un conjunto de átomos en una región específica del sistema.

    Parameters
    ----------
    file : str
        Ruta del archivo que contiene la simulación.
    atom : str
        Nombre del átomo a seleccionar.
    z_lo : float
        Valor inferior del rango en el eje Z.
    z_hi : float
        Valor superior del rango en el eje Z.

    Returns
    -------
    np.ndarray
        Array con los valores de MSD para cada tiempo.
    np.ndarray
        Array con los valores de tiempo de lag correspondientes a cada MSD.

    """
    # Carga la trayectoria y selecciona los átomos en la región específica
    u = mda.Universe(file)
    water = u.select_atoms("prop {} < z and prop z < {}".format(z_lo, z_hi))

    # Calcula el MSD utilizando la clase EinsteinMSD de MDAnalysis
    MSD = msda.EinsteinMSD(water, select=atom, msd_type='xyz', fft=True)
    MSD.run()

    # Obtiene los resultados de MSD y los tiempos de lag
    msd = MSD.results.timeseries
    nframes = MSD.n_frames
    timestep = 1  # Esto debe ser el tiempo real entre frames
    lagtimes = np.arange(nframes) * timestep

    return msd, lagtimes

def plot_msd(index_start, index_end, msd, lagtimes):
    """
    Grafica el MSD y calcula el coeficiente de difusión.

    Parameters
    ----------
    index_start : int
        Índice de inicio para el cálculo del MSD y la regresión lineal.
    index_end : int
        Índice de fin para el cálculo del MSD y la regresión lineal.
    msd : np.ndarray
        Array con los valores de MSD para cada tiempo.
    lagtimes : np.ndarray
        Array con los valores de tiempo de lag correspondientes a cada MSD.

    Returns
    -------
    float
        Coeficiente de difusión calculado a partir de la pendiente de la regresión lineal.

    """
    # Realiza la regresión lineal del MSD en el intervalo especificado
    slope, intercept, r_value, p_value, std_err = stats.linregress(lagtimes[index_start:index_end], msd[index_start:index_end])

    # Grafica el MSD con la regresión lineal
    ax = sns.regplot(x=lagtimes[index_start:index_end], y=msd[index_start:index_end], line_kws={'label': "y={0:.5f}x+{1:.5f}".format(slope, intercept)})
    ax.set(xlabel="MSD (A²)", ylabel="Time (ps)")
    ax.legend()
    plt.show()

    # Calcula el coeficiente de difusión
    diffusion = slope / 6

    return diffusion

def msd_z_grid(atom, file, image_name):
    """
    Calcula y grafica el coeficiente de difusión del agua en diferentes regiones a lo largo del eje Z.

    Parameters
    ----------
    atom : str
        Nombre del átomo a seleccionar.
    file : str
        Ruta del archivo que contiene la simulación.
    image_name : str
        Nombre del archivo de imagen para guardar la visualización.

    """
    # Crear un DataFrame vacío para almacenar los resultados
    df = pd.DataFrame(columns=["z", "diffusion"])

    # Crear un array de numpy con los valores de z_hi directamente
    z = np.arange(0, 100, 10)  # Modificar según las regiones deseadas

    # Crear un array de numpy para almacenar los resultados de msd_z y plot_msd
    diffusion = np.zeros(len(z))

    # Iterar sobre los valores de z y calcular msd y diffusion
    for i, z_hi in enumerate(z):
        z_lo = z_hi - 5
        m = msd_z(file, atom, z_lo, z_hi)
        diffusion[i] = plot_msd(10, 60, m[0], m[1])

    # Crear un DataFrame de pandas con los resultados
    df = pd.DataFrame({'z': z, 'diffusion': diffusion * 10})

    # Graficar los resultados
    dosd_plot(df.diffusion, df.z, image_name)

def dosd_plot(x, y, image_name):
    """
    Grafica la densidad del coeficiente de difusión en función de la posición Z.

    Parameters
    ----------
    x : np.ndarray
        Array con los valores de coeficiente de difusión.
    y : np.ndarray
        Array con los valores de posición Z.
    image_name : str
        Nombre del archivo de imagen para guardar la visualización.

    """
    sns.set_style("white")
    sns.kdeplot(x=x, y=y, cmap="BuPu", shade=False, bw_adjust=0.90)

    # Configuración del eje X
    plt.rcParams['font.family'] = 'Arial'
    plt.xlabel('Difusión / x10$^{-9}$ m$^2$.s$^{-1}$', fontsize=20)
    plt.xticks(fontsize=13)
    plt.xlim(0, 2.5)
    plt.xticks(np.arange(0.0, 3.1, 0.5), fontsize=15)
    plt.xticks(np.arange(0.25, 3.1, 0.25), minor=True)
    plt.tick_params(axis='x', which='both', bottom=True, top=True)
    plt.gca().xaxis.set_tick_params(which='both', width=1, direction='in', length=7, pad=5)

    # Configuración del eje Y
    plt.ylabel('Posición Z (Å)', fontsize=20)
    plt.ylim(0, 100)
    plt.yticks(np.arange(25, 125, 25), fontsize=15)
    plt.yticks(np.arange(12.5, 112.5, 12.5), minor=True)
    plt.tick_params(axis='y', which='both', left=True, right=True)
    plt.gca().yaxis.set_tick_params(which='both', width=1, direction='in', length=7, pad=5)

    plt.rcParams['axes.linewidth'] = 1

    plt.tight_layout()
    plt.savefig(image_name, dpi=200)
    plt.show()


        
    
    
    

