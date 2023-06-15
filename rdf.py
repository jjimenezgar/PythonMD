import argparse
import numpy as np
import matplotlib.pyplot as plt
import MDAnalysis as mda
from MDAnalysis.analysis import rdf

def read_input_file(input_file):
    """
    Lee un archivo de entrada y crea una instancia de Universe.

    Args:
        input_file (str): Ruta del archivo de entrada.

    Returns:
        Universe: Instancia de Universe creada a partir del archivo de entrada.

    Raises:
        ValueError: Si no se puede leer el archivo de entrada.
    """
    try:
        u = mda.Universe(input_file)
    except Exception as e:
        raise ValueError(f'Error: Could not read input file {input_file}: {e}')
    return u

def create_large_cell(u):
    """
    Crea una celda grande duplicando las dimensiones de la celda original.

    Args:
        u (Universe): Instancia de Universe.

    Modifica:
        Modifica las dimensiones de la instancia de Universe para crear una celda grande.
    """
    d = np.max(u.atoms.positions) - np.min(u.atoms.positions)
    u.dimensions = [2*d, 2*d, 2*d, 90, 90, 90]

def calculate_rdf(u, atom_sel1, atom_sel2):
    """
    Calcula la función de distribución radial (RDF) para dos selecciones de átomos.

    Args:
        u (Universe): Instancia de Universe.
        atom_sel1 (str): Selección de átomos 1.
        atom_sel2 (str): Selección de átomos 2.

    Returns:
        ndarray, ndarray: Arreglos de NumPy que contienen los valores del RDF y los bins correspondientes.

    Nota:
        Las selecciones de átomos deben seguir la sintaxis utilizada por MDAnalysis.

    Referencia:
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
    Guarda los datos del RDF en un archivo de salida.

    Args:
        output_file (str): Ruta del archivo de salida.
        bins (ndarray): Arreglo de NumPy que contiene los bins del RDF.
        d_rdf (ndarray): Arreglo de NumPy que contiene los valores del RDF.

    Nota:
        Los datos se guardan en formato CSV con los bins y los valores del RDF como columnas.
    """
    data = np.column_stack((bins, d_rdf))
    np.savetxt(output_file, data, delimiter=',')

def plot_rdf(bins, d_rdf):
    """
    Grafica el RDF.

    Args:
        bins (ndarray): Arreglo de NumPy que contiene los bins del RDF.
        d_rdf (ndarray): Arreglo de NumPy que contiene los valores del RDF.
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
    Procesa el cálculo del RDF.

    Args:
        input_file (str): Ruta del archivo de entrada.
        output_file (str): Ruta del archivo de salida.
        atom_sel1 (str): Selección de átomos 1.
        atom_sel2 (str): Selección de átomos 2.
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
