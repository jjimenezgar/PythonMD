#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This script calculates the coordination number of a set of oxygen atoms
in a molecular dynamics simulation using the MDAnalysis library.
"""

import MDAnalysis as mda
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


def cn_atom(file, r_cut_lo, r_cut_hi):
    """
    Calculates the coordination number of a set of oxygen atoms.

    Parameters
    ----------
    file : str
        Path to the simulation file.
    r_cut_lo : float
        Lower cutoff radius to search for oxygen atoms.
    r_cut_hi : float
        Upper cutoff radius to search for oxygen atoms.

    Returns
    -------
    float
        The average coordination number of the set of oxygen atoms.

    """
    # Load the trajectory and topology
    u = mda.Universe(file)

    # Select only the oxygen atoms
    oxygen = u.select_atoms("prop 30 < z and prop z < 60 and name O1 and around {} name O1 and not around {} name O1".format(r_cut_hi, r_cut_lo))
    
    coord = [len(u.select_atoms("name O1 and around 3.25 index {}".format(atom.index))) 
         for ts in u.trajectory 
         if ts.frame == 0 
         for atom in oxygen]
         
    return np.mean(coord)


def graf_cn(X, Y):
    """
    Plots the coordination number as a function of radius.

    Parameters
    ----------
    X : array-like
        Array with the radii.
    Y : array-like
        Array with the corresponding coordination numbers.

    Returns
    -------
    None.

    """
    plt.plot(X, Y, color='black', linestyle='dashed', linewidth=2, label="12")

    plt.title('Coordination number/water', fontsize=18)

    # Add labels to the axes
    plt.xlabel('Radius (Å)', fontsize=16)
    plt.ylabel('Coordination Number (CN)', fontsize=16)

    plt.legend(fontsize=12)

    # Show the plot
    plt.show()


def processing(olig):
    """
    Processes the simulation to calculate the average coordination number
    of a set of oxygen atoms.

    Parameters
    ----------
    olig : int
        Oligomer number to process.

    Returns
    -------
    None.

    """
    # Define the simulation file path
    file = "~/{}_cn.xyz".format(olig)

    # Create an empty DataFrame to store the results
    df = pd.DataFrame(columns=["R", "CN"])

    # Define parameters for coordination number calculation
    radius_p = 20
    step = 5
    final = radius_p + step

    # Iterates over the r_cut values to calculate the coordination number
    for i in np.arange(step, final, step):
        # Calculate the value of r_cut
        r_cut_lo = (final) - i
        r_cut_hi = radius_p - i

        # Calculate the coordination number
        CN = cn_atom(file, r_cut_hi, r_cut_lo)

        # Add a row to the DataFrame with the r_cut values and the current coordination number
        df = df.append({"R": r_cut_lo, "CN": CN}, ignore_index=True)

    # Plot the results
    graf_cn(df.R, df.CN)


if __name__ == '__main__':
    for olig in [4, 8, 12, 24, 36]:
        processing(olig)
