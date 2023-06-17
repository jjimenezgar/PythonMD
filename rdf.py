import argparse
import numpy as np
import matplotlib.pyplot as plt
import MDAnalysis as mda
from MDAnalysis.analysis import rdf

def read_input_file(input_file):
    """
    Read an input file and create a Universe instance.

    Args:
        input_file (str): Path to the input file.

    Returns:
        Universe: Universe instance created from the input file.

    Raises:
        ValueError: If the input file cannot be read.
    """
    try:
        u = mda.Universe(input_file)
    except Exception as e:
        raise ValueError(f'Error: Could not read input file {input_file}: {e}')
    return u

def create_large_cell(u):
    """
    Create a large cell by doubling the dimensions of the original cell.

    Args:
        u (Universe): Universe instance.

    Modifies:
        Modifies the dimensions of the Universe instance to create a large cell.
    """
    d = np.max(u.atoms.positions) - np.min(u.atoms.positions)
    u.dimensions = [2*d, 2*d, 2*d, 90, 90, 90]

def calculate_rdf(u, atom_sel1, atom_sel2):
    """
    Calculate the radial distribution function (RDF) for two atom selections.

    Args:
        u (Universe): Universe instance.
        atom_sel1 (str): Atom selection 1.
        atom_sel2 (str): Atom selection 2.

    Returns:
        ndarray, ndarray: NumPy arrays containing the RDF values and the corresponding bins.

    Note:
        The atom selections must follow the syntax used by MDAnalysis.

    Reference:
        https://docs.mdanalysis.org/stable/documentation_pages/selections.html
    """
    atom_1 = u.select_atoms(atom_sel1)
    atom_2 = u.select_atoms(atom_sel2)
    
    atoms_rdf = rdf.InterRDF(atom_1, atom_2)
    atoms_rdf.run()
    
    density = u.atoms.n_atoms / mda.lib.mdamath.box_volume(u.dimensions)
    d_rdf = density * atoms_rdf.rdf
    bins = atoms_rdf.bins
    
    return d_rdf, bins

def save_data(output_file, bins, d_rdf):
    """
    Save RDF data to an output file.

    Args:
        output_file (str): Path to the output file.
        bins (ndarray): NumPy array containing the RDF bins.
        d_rdf (ndarray): NumPy array containing the RDF values.

    Note:
        The data is saved in CSV format with the bins and RDF values as columns.
    """
    data = np.column_stack((bins, d_rdf))
    np.savetxt(output_file, data, delimiter=',')

def plot_rdf(bins, d_rdf):
    """
    Plot the RDF.

    Args:
        bins (ndarray): NumPy array containing the RDF bins.
        d_rdf (ndarray): NumPy array containing the RDF values.
    """
    fig, ax = plt.subplots()
    ax.set_xlim(1, 10)
    ax.set_ylim(0, 0.8)
    
    ax.plot(bins, d_rdf)
    plt.title('RDF', fontsize=24)
    plt.xlabel('r (Å)', fontsize=22)
    plt.xticks(fontsize=12)
    plt.ylabel('g(r)', fontsize=22)
    plt.yticks(fontsize=12)
    plt.show()

def process_rdf(input_file, output_file, atom_sel1, atom_sel2):
    """
    Process the RDF calculation.

    Args:
        input_file (str): Path to the input file.
        output_file (str): Path to the output file.
        atom_sel1 (str): Atom selection 1.
        atom_sel2 (str): Atom selection 2.
    """
    # Parse input file
    u = read_input_file(input_file)

    # Create large cell
    create_large_cell(u)

    # Calculate RDF
    d_rdf, bins = calculate_rdf(u, atom_sel1, atom_sel2)

    # Save data
    save_data(output_file, bins, d_rdf)

    # Plot RDF
    plot_rdf(bins, d_rdf)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Calculate RDF')
    parser.add_argument('input_file', type=str, help='Input file path')
    parser.add_argument('output_file', type=str, help='Output file path')
    parser.add_argument('atom_sel1', type=str, help='Atom selection 1')
    parser.add_argument('atom_sel2', type=str, help='Atom selection 2')
    args = parser.parse_args()

    process_rdf(args.input_file, args.output_file, args.atom_sel1, args.atom_sel2)
