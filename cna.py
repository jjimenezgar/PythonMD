import numpy as np
from ovito.io import import_file
from ovito.modifiers import CommonNeighborAnalysisModifier, SliceModifier, SelectTypeModifier, DeleteSelectedModifier, ColorCodingModifier
from ovito.vis import Viewport, OSPRayRenderer

def calculate_atomic_strain(input_file):
    """
    Calcula la deformación atómica y la estructura del sistema a partir de un archivo XYZ y guarda una imagen de la simulación.

    Args:
        input_file (str): Ruta del archivo XYZ de entrada.

    Returns:
        None
    """
    # Importar archivo
    pipeline = import_file(input_file)
    pipeline.add_to_scene()
    data = pipeline.compute()

    # Aplicar el modificador SelectTypeModifier para seleccionar todas las partículas del tipo 'Ox'
    pipeline.modifiers.append(SelectTypeModifier(
        operate_on="particles",
        property="Particle Type",
        types={'O1', 'O2', 'O3', 'O4', 'H1', 'H2', 'C', 'Cgr', 'F', 'S'}
    ))
    # Agregar el modificador DeleteSelectedModifier para eliminar las partículas seleccionadas
    pipeline.modifiers.append(DeleteSelectedModifier())

    # Agregar el modificador CommonNeighborAnalysisModifier para determinar el tipo estructural de cada átomo
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

    print("Fracción de FCC:", pipeline.compute().attributes['fcc_fraction'])
    print("Fracción de HCP:", pipeline.compute().attributes['hcp_fraction'])
    print("Fracción de otros:", pipeline.compute().attributes['other_fraction'])

    # Agregar el modificador SliceModifier para cortar la celda
    pipeline.modifiers.append(SliceModifier(
        distance=220.0  # Este valor está relacionado con la posición donde hacer el corte en la celda
    ))

    # Agregar el modificador ColorCodingModifier para codificar el color según el tipo de estructura
    pipeline.modifiers.append(ColorCodingModifier(
        operate_on="particles",
        property='Structure Type',
        gradient=ColorCodingModifier.Rainbow()
    ))

    data.cell.vis.enabled = False  # Desactivar la visualización de la celda

    # Renderizar y guardar una imagen de la simulación
    vp = Viewport(type=Viewport.Type.Right)
    vp.zoom_all()
    vp.render_image(filename='simulation.png',
                    size=(1000, 1000),
                    renderer=OSPRayRenderer())

    print("La simulación se ha guardado como 'simulation.png'.")

# Ejemplo de uso
input_file = "/path/to/input_file.xyz"
calculate_atomic_strain(input_file)
