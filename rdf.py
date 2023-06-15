import argparse
import numpy as np
import matplotlib.pyplot as plt
import MDAnalysis as mda
from MDAnalysis.analysis import rdf

def read_input_file(input_file):
    try:
        u = mda.Universe(input_file)
    except Exception as e:
        raise ValueError(f'Error: Could not read input file {input_file}: {e}')
    return u

def create_large_cell(u):
    d = np.max(u.atoms.positions) - np.min(u.atoms.positions)
    u.dimensions = [2*d, 2*d, 2*d, 90, 90, 90]

def calculate_rdf(u, atom_sel1, atom_sel2):
    atom_1 = u.select_atoms(atom_sel1)
    atom_2 = u.select_atoms(atom_sel2)
    
    atoms_rdf = rdf.InterRDF(atom_1, atom_2)
    atoms_rdf.run()
    
    density = u.atoms.n_atoms / mda.lib.mdamath.box_volume(u.dimensions)
    d_rdf = density * atoms_rdf.rdf
    bins = atoms_rdf.bins
    
    return d_rdf, bins

def save_data(output_file, bins, d_rdf):
    data = np.column_stack((bins, d_rdf))
    np.savetxt(output_file, data, delimiter=',')

def plot_rdf(bins, d_rdf):
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
