import numpy as np
from ovito.io import import_file
from ovito.modifiers import AtomicStrainModifier, SelectTypeModifier, DeleteSelectedModifier, ColorCodingModifier
from ovito.vis import Viewport, OSPRayRenderer

def calculate_atomic_strain(input_file):
    """
    Calcula la deformación atómica de un sistema y guarda una imagen de la simulación.

    Args:
        input_file (str): Ruta del archivo XYZ de entrada.

    Returns:
        None
    """
    # Importar archivo
    pipeline = import_file(input_file)
    pipeline.add_to_scene()
    data = pipeline.compute()

    # Seleccionar todas las partículas que no queremos estudiar
    pipeline.modifiers.append(SelectTypeModifier(
        operate_on="particles",
        property="Particle Type",
        types={'Atoms to delete'}
    ))
    # Eliminar las partículas seleccionadas
    pipeline.modifiers.append(DeleteSelectedModifier())

    # Calcular la deformación atómica
    pipeline.modifiers.append(AtomicStrainModifier(
        cutoff=3.0
    ))

    # Codificar el color según la deformación de corte
    pipeline.modifiers.append(ColorCodingModifier(
        operate_on="particles",
        property='Shear Strain',
        gradient=ColorCodingModifier.Rainbow()
    ))

    data.cell.vis.enabled = False  # Desactivar la visualización de la celda

    # Renderizar y guardar una imagen de la simulación
    vp = Viewport(type=Viewport.Type.Front)
    vp.zoom_all()
    vp.render_image(filename='simulation.png',
                    size=(1000, 1000),
                    renderer=OSPRayRenderer())

    print("La simulación se ha guardado como 'simulation.png'.")

# Ejemplo de uso
input_file = "/path/to/input_file.xyz"
calculate_atomic_strain(input_file)
