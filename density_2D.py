import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis import density
import matplotlib.pyplot as plt

# Input file
file = "input.xyz"

# Create Universe object
u = mda.Universe(file)

def density_2D(u):
    """
    Calculates and visualizes the 2D average density of water in a simulation box.

    Args:
        u (Universe): MDAnalysis Universe object.

    Returns:
        None
    """

    # Define the maximum Z coordinate of the box
    z_hi = "Enter maximum box Z value"

    # Select oxygen atoms in the desired region (above z_hi)
    ow = u.select_atoms("name O1 and prop {} > z".format(z_hi))

    # Perform density analysis
    dens = density.DensityAnalysis(ow, delta=3.5, padding=0)
    dens.run()

    # Get the density grid and convert it to SPC/E units
    grid = dens.results.density.grid
    dens.results.density.convert_density('SPC')

    # Calculate the average density along the Z axis (top view)
    avg = grid.mean(axis=-1)

    # Create the figure and axes
    fig, ax = plt.subplots(figsize=(10, 8))

    # Create the heat map of the average density
    im = ax.imshow(avg, interpolation="bicubic", cmap='jet', vmin=0, vmax=1)

    # Add a colorbar
    cbar = plt.colorbar(im)
    cbar.set_label('Mean density of water over SPC/E literature value')

    # Axis labels
    plt.xlabel('X-axis ($\AA$)')
    plt.ylabel('Y-axis ($\AA$)')

# Call the function to perform the analysis
density_2D(u)
