import numpy as np
from ovito.io import import_file, export_file
from ovito.modifiers import CoordinationAnalysisModifier, SelectTypeModifier, DeleteSelectedModifier
from ovito.vis import Viewport, OSPRayRenderer


def perform_cluster_analysis(input_file, particle_types, cutoff=3.0):
    """
    Performs coordination analysis on a given input file for specific particle types.

    :param input_file: Path to the input file containing particle data.
    :param particle_types: Set of particle types to include in the analysis.
    :param cutoff: Distance cutoff for coordination determination (default: 3.0).
    """
    # Import file
    pipeline = import_file(input_file)
    pipeline.add_to_scene()

    # Compute pipeline data
    data = pipeline.compute()

    # Apply SelectTypeModifier to select particles of specific types
    pipeline.modifiers.append(SelectTypeModifier(
        operate_on="particles",
        property="Particle Type",
        types=particle_types
    ))

    # Apply DeleteSelectedModifier to remove selected particles
    pipeline.modifiers.append(DeleteSelectedModifier())

    # Apply coordinate analysis
    pipeline.modifiers.append(CoordinationAnalysisModifier(
        cutoff=cutoff))

    # Export results of the coordinate analysis
    export_file(pipeline, "/outputfile.xyz", "xyz",
    columns = [ "Coordination","Position.X", "Position.Y", "Position.Z"])


    # Disable cell visualization
    data.cell.vis.enabled = False

    # Configure viewport and render image
    vp = Viewport(type=Viewport.Type.Right)
    vp.zoom_all()
    vp.render_image(filename='simulation.png', size=(1000, 1000), renderer=OSPRayRenderer())


if __name__ == "__main__":
    input_file = "/home/jjimenez/Documentos/ApendiceA/cn.xyz"
    particle_types = {'C Cgr O1 O2 O3 O4 S H1 H2 F'}  # Set of particle types to exclude in the analysis
    perform_cluster_analysis(input_file, particle_types)