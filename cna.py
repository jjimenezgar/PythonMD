import numpy as np
from ovito.io import import_file
from ovito.modifiers import CommonNeighborAnalysisModifier, SliceModifier, SelectTypeModifier, DeleteSelectedModifier, ColorCodingModifier
from ovito.vis import Viewport, OSPRayRenderer

def calculate_atomic_strain(input_file):
    """
    Calculates atomic strain and structure of the system from an XYZ file and saves an image of the simulation.

    Args:
        input_file (str): Path to the input XYZ file.

    Returns:
        None
    """
    # Import the file
    pipeline = import_file(input_file)
    pipeline.add_to_scene()
    data = pipeline.compute()

    # Apply SelectTypeModifier to select all particles of type 'Ox'
    pipeline.modifiers.append(SelectTypeModifier(
        operate_on="particles",
        property="Particle Type",
        types={'O1', 'O2', 'O3', 'O4', 'H1', 'H2', 'C', 'Cgr', 'F', 'S'}
    ))
    # Add DeleteSelectedModifier to delete the selected particles
    pipeline.modifiers.append(DeleteSelectedModifier())

    # Add CommonNeighborAnalysisModifier to determine the structural type of each atom
    pipeline.modifiers.append(CommonNeighborAnalysisModifier())

    def compute_fraction(frame, data):
        n_fcc = data.attributes['CommonNeighborAnalysis.counts.FCC']
        n_hcp = data.attributes['CommonNeighborAnalysis.counts.HCP']
        n_other = data.attributes['CommonNeighborAnalysis.counts.OTHER']
        total_count = data.particles.count

        data.attributes['fcc_fraction'] = n_fcc / total_count
        data.attributes['hcp_fraction'] = n_hcp / total_count
        data.attributes['other_fraction'] = n_other / total_count

    pipeline.modifiers.append(compute_fraction)

    print("FCC Fraction:", pipeline.compute().attributes['fcc_fraction'])
    print("HCP Fraction:", pipeline.compute().attributes['hcp_fraction'])
    print("Other Fraction:", pipeline.compute().attributes['other_fraction'])

    # Add SliceModifier to cut the cell
    pipeline.modifiers.append(SliceModifier(
        distance=220.0  # This value is related to the position to make the cut in the cell
    ))

    # Add ColorCodingModifier to color code the particles based on the structure type
    pipeline.modifiers.append(ColorCodingModifier(
        operate_on="particles",
        property='Structure Type',
        gradient=ColorCodingModifier.Rainbow()
    ))

    data.cell.vis.enabled = False  # Disable cell visualization

    # Render and save an image of the simulation
    vp = Viewport(type=Viewport.Type.Right)
    vp.zoom_all()
    vp.render_image(filename='simulation.png',
                    size=(1000, 1000),
                    renderer=OSPRayRenderer())

    print("The simulation has been saved as 'simulation.png'.")

if __name__ == '__main__':
    # Example usage
    input_file = input("Enter the path to the file: ")
    calculate_atomic_strain(input_file)
