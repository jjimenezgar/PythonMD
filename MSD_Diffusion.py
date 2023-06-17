import MDAnalysis.analysis.msd as msda
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
import seaborn as sns
import pandas as pd

def msd_z(file, atom, z_lo, z_hi):
    """
    Calculates the Mean Squared Displacement (MSD) of a set of atoms in a specific region of the system.

    Parameters
    ----------
    file : str
        Path to the simulation file.
    atom : str
        Name of the atom to select.
    z_lo : float
        Lower value of the range on the Z-axis.
    z_hi : float
        Upper value of the range on the Z-axis.

    Returns
    -------
    np.ndarray
        Array with the MSD values for each time.
    np.ndarray
        Array with the corresponding lag times for each MSD.

    """
    # Load the trajectory and select atoms in the specific region
    u = mda.Universe(file)
    water = u.select_atoms("prop {} < z and prop z < {}".format(z_lo, z_hi))

    # Calculate MSD using the EinsteinMSD class from MDAnalysis
    MSD = msda.EinsteinMSD(water, select=atom, msd_type='xyz', fft=True)
    MSD.run()

    # Get the MSD results and lag times
    msd = MSD.results.timeseries
    nframes = MSD.n_frames
    timestep = 1  # This should be the actual time between frames
    lagtimes = np.arange(nframes) * timestep

    return msd, lagtimes

def plot_msd(index_start, index_end, msd, lagtimes):
    """
    Plots the MSD and calculates the diffusion coefficient.

    Parameters
    ----------
    index_start : int
        Start index for calculating MSD and performing linear regression.
    index_end : int
        End index for calculating MSD and performing linear regression.
    msd : np.ndarray
        Array with the MSD values for each time.
    lagtimes : np.ndarray
        Array with the corresponding lag times for each MSD.

    Returns
    -------
    float
        Diffusion coefficient calculated from the slope of the linear regression.

    """
    # Perform linear regression of MSD within the specified interval
    slope, intercept, r_value, p_value, std_err = stats.linregress(lagtimes[index_start:index_end], msd[index_start:index_end])

    # Plot MSD with linear regression
    ax = sns.regplot(x=lagtimes[index_start:index_end], y=msd[index_start:index_end], line_kws={'label': "y={0:.5f}x+{1:.5f}".format(slope, intercept)})
    ax.set(xlabel="MSD (A²)", ylabel="Time (ps)")
    ax.legend()
    plt.show()

    # Calculate the diffusion coefficient
    diffusion = slope / 6

    return diffusion

def msd_z_grid(atom, file, image_name):
    """
    Calculates and plots the diffusion coefficient of water in different regions along the Z-axis.

    Parameters
    ----------
    atom : str
        Name of the atom to select.
    file : str
        Path to the simulation file.
    image_name : str
        Name of the image file to save the visualization.

    """
    # Create an empty DataFrame to store the results
    df = pd.DataFrame(columns=["z", "diffusion"])

    # Create a numpy array with the z_hi values directly
    z = np.arange(0, 100, 10)  # Modify according to desired regions

    # Create a numpy array to store the results of msd_z and plot_msd
    diffusion = np.zeros(len(z))

    # Iterate over z values and calculate msd and diffusion
    for i, z_hi in enumerate(z):
        z_lo = z_hi - 5
        m = msd_z(file, atom, z_lo, z_hi)
        diffusion[i] = plot_msd(10, 60, m[0], m[1])

    # Create a pandas DataFrame with the results
    df = pd.DataFrame({'z': z, 'diffusion': diffusion * 10})

    # Plot the results
    dosd_plot(df.diffusion, df.z, image_name)

def dosd_plot(x, y, image_name):
    """
    Plots the density of diffusion coefficient as a function of the Z position.

    Parameters
    ----------
    x : np.ndarray
        Array with the diffusion coefficient values.
    y : np.ndarray
        Array with the Z position values.
    image_name : str
        Name of the image file to save the visualization.

    """
    sns.set_style("white")
    sns.kdeplot(x=x, y=y, cmap="BuPu", shade=False, bw_adjust=0.90)

    # X-axis configuration
    plt.rcParams['font.family'] = 'Arial'
    plt.xlabel('Diffusion / x10$^{-9}$ m$^2$.s$^{-1}$', fontsize=20)
    plt.xticks(fontsize=13)
    plt.xlim(0, 2.5)
    plt.xticks(np.arange(0.0, 3.1, 0.5), fontsize=15)
    plt.xticks(np.arange(0.25, 3.1, 0.25), minor=True)
    plt.tick_params(axis='x', which='both', bottom=True, top=True)
    plt.gca().xaxis.set_tick_params(which='both', width=1, direction='in', length=7, pad=5)

    # Y-axis configuration
    plt.ylabel('Z Position (Å)', fontsize=20)
    plt.ylim(0, 100)
    plt.yticks(np.arange(25, 125, 25), fontsize=15)
    plt.yticks(np.arange(12.5, 112.5, 12.5), minor=True)
    plt.tick_params(axis='y', which='both', left=True, right=True)
    plt.gca().yaxis.set_tick_params(which='both', width=1, direction='in', length=7, pad=5)

    plt.rcParams['axes.linewidth'] = 1

    plt.tight_layout()
    plt.savefig(image_name, dpi=200)
    plt.show()
