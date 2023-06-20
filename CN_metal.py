import sys
import numpy as np
from ovito.io import import_file, export_file
from ovito.modifiers import CoordinationAnalysisModifier, SelectTypeModifier, DeleteSelectedModifier
from ovito.vis import Viewport, OSPRayRenderer

def perform_cluster_analysis(input_file, excluded_particles, cutoff=3.0):
    """
    Performs coordination analysis on a given input file for specific particle types.

    :param input_file: Path to the input file containing particle data.
    :param excluded_particles: Set of particle types to exclude from the analysis.
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
        types=~excluded_particles
    ))

    # Apply DeleteSelectedModifier to remove selected particles
    pipeline.modifiers.append(DeleteSelectedModifier())

    # Apply coordinate analysis
    pipeline.modifiers.append(CoordinationAnalysisModifier(cutoff=cutoff))

    # Export results of the coordinate analysis
    export_file(pipeline, "/metal.xyz", "xyz",
                columns=["Coordination", "Position.X", "Position.Y", "Position.Z"])

    # Disable cell visualization
    data.cell.vis.enabled = False

    # Configure viewport and render image
    vp = Viewport(type=Viewport.Type.Right)
    vp.zoom_all()
    vp.render_image(filename='simulation.png', size=(1000, 1000), renderer=OSPRayRenderer())

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python script.py input_file excluded_particles")
        sys.exit(1)

    input_file = sys.argv[1]
    excluded_particles = set(sys.argv[2].split())

    perform_cluster_analysis(input_file, excluded_particles)
