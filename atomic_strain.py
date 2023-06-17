import numpy as np
from ovito.io import import_file
from ovito.modifiers import AtomicStrainModifier, SelectTypeModifier, DeleteSelectedModifier, ColorCodingModifier
from ovito.vis import Viewport, OSPRayRenderer

def calculate_atomic_strain(input_file):
    """
    Calculates the atomic strain of a system and saves an image of the simulation.

    Args:
        input_file (str): Path to the input XYZ file.

    Returns:
        None
    """
    # Import file
    pipeline = import_file(input_file)
    pipeline.add_to_scene()
    data = pipeline.compute()

    # Select all particles we don't want to study
    pipeline.modifiers.append(SelectTypeModifier(
        operate_on="particles",
        property="Particle Type",
        types={'Atoms to delete'}
    ))
    # Delete selected particles
    pipeline.modifiers.append(DeleteSelectedModifier())

    # Calculate atomic strain
    pipeline.modifiers.append(AtomicStrainModifier(
        cutoff=3.0
    ))

    # Color code based on shear strain
    pipeline.modifiers.append(ColorCodingModifier(
        operate_on="particles",
        property='Shear Strain',
        gradient=ColorCodingModifier.Rainbow()
    ))

    data.cell.vis.enabled = False  # Disable cell visualization

    # Render and save an image of the simulation
    vp = Viewport(type=Viewport.Type.Front)
    vp.zoom_all()
    vp.render_image(filename='simulation.png',
                    size=(1000, 1000),
                    renderer=OSPRayRenderer())

    print("The simulation has been saved as 'simulation.png'.")

if __name__ == "__main__":
    # Example usage
    input_file = "/path/to/input_file.xyz"
    calculate_atomic_strain(input_file)
