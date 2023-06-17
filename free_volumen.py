import numpy as np
from ovito.io import import_file
from ovito.modifiers import ConstructSurfaceModifier, SelectTypeModifier, DeleteSelectedModifier, AffineTransformationModifier

def calculate_free_volume(input_file, lz):
    """
    Calculates the free volume of a molecular system from an input file in XYZ format.

    Parameters:
    - input_file (str): Path to the input file in XYZ format containing the molecular structure.
    - lz (float): Value of the z coordinate of the cell in angstroms. It is used to modify the system's cell.

    Returns:
    - free_volume (float): Free volume of the system in cubic meters.

    Requires the OVITO library: https://www.ovito.org/

    The function performs the following steps:
    1. Imports the specified XYZ file.
    2. Modifies the cell vectors of the system to adjust the z coordinate value to 'lz'.
    3. Applies transformations to the cell using the AffineTransformationModifier in OVITO.
    4. Selects particles of type 'Ox' using the SelectTypeModifier.
    5. Deletes the selected particles using the DeleteSelectedModifier.
    6. Constructs a surface based on the AlphaShape method using the ConstructSurfaceModifier.
    7. Calculates the free volume by subtracting the volume of the empty cell from the volume of the cell occupied by the system.
    8. Converts the volume from cubic angstroms to cubic meters.

    Example usage:
    input_file = "/path/to/file.xyz"
    lz = 60.0
    free_volume, total_volume = calculate_free_volume(input_file, lz)
    print(f"Free volume: {free_volume} m³")
    """

    # Import file
    pipeline = import_file(input_file)
    data = pipeline.compute()

    lx = data.cell[0, 0]
    ly = data.cell[1, 1]

    # Modify the cell vectors
    data.cell[2, 2] = lz
    data.cell.pbc = (True, True, True)

    # Add the AffineTransformationModifier to apply absolute transformations to the cell
    pipeline.modifiers.append(AffineTransformationModifier(
        relative_mode=False,
        target_cell=pipeline.compute(0).cell[...]
    ))

    # Apply the SelectTypeModifier to select all particles of type 'Ox'
    pipeline.modifiers.append(SelectTypeModifier(
        operate_on="particles",
        property="Particle Type",
        types={'Ox', 'Cgr'}
    ))

    # Add the AffineTransformationModifier to apply absolute transformations to the cell again
    pipeline.modifiers.append(AffineTransformationModifier(
        relative_mode=False,
        target_cell=data.cell[...]
    ))

    # Add the DeleteSelectedModifier to delete the selected particles
    pipeline.modifiers.append(DeleteSelectedModifier())

    # Add the ConstructSurfaceModifier to construct a surface based on the AlphaShape method
    pipeline.modifiers.append(ConstructSurfaceModifier(
        method=ConstructSurfaceModifier.Method.AlphaShape,
        radius=3.0,
        identify_regions=True
    ))

    data = pipeline.compute()

    empty_volume = data.attributes['ConstructSurfaceMesh.empty_volume']

    total_volume = (lx * ly * lz) * 1e-30  # Convert from A^3 to m^3
    free_volume = empty_volume * 1e-30  # Convert from A^3 to m^3

    return free_volume, total_volume
