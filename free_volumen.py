import numpy as np
from ovito.io import import_file
from ovito.modifiers import ConstructSurfaceModifier, SelectTypeModifier, DeleteSelectedModifier, AffineTransformationModifier

def calculate_free_volume(input_file, lz):
    """
    Calcula el volumen libre de un sistema molecular a partir de un archivo de entrada en formato XYZ.

    Parámetros:
    - input_file (str): Ruta al archivo de entrada en formato XYZ que contiene la estructura molecular.
    - lz (float): Valor de la coordenada z de la celda en angstroms. Se utiliza para modificar la celda del sistema.

    Retorna:
    - free_volume (float): Volumen libre del sistema en metros cúbicos.

    Requiere la biblioteca OVITO: https://www.ovito.org/

    La función realiza los siguientes pasos:
    1. Importa el archivo XYZ especificado.
    2. Modifica los vectores de celda del sistema para ajustar el valor de la coordenada z a 'lz'.
    3. Aplica transformaciones a la celda utilizando el modificador AffineTransformationModifier de OVITO.
    4. Selecciona las partículas del tipo 'Ox' utilizando el modificador SelectTypeModifier.
    5. Elimina las partículas seleccionadas utilizando el modificador DeleteSelectedModifier.
    6. Construye una superficie basada en el método AlphaShape utilizando el modificador ConstructSurfaceModifier.
    7. Calcula el volumen libre restando el volumen de la celda vacía al volumen de la celda ocupada por el sistema.
    8. Convierte el volumen de angstrom cúbico a metros cúbicos.

    Ejemplo de uso:
    input_file = "/ruta/al/archivo.xyz"
    lz = 60.0
    free_volume, total_volume = calculate_free_volume(input_file, lz)
    print(f"Volumen libre: {free_volume} m³")
    """

    # Importar archivo
    pipeline = import_file(input_file)
    data = pipeline.compute()

    lx = data.cell[0, 0]
    ly = data.cell[1, 1]

    # Modificar los vectores de celda
    data.cell[2, 2] = lz
    data.cell.pbc = (True, True, True)

    # Agregar el modificador AffineTransformationModifier para aplicar transformaciones absolutas a la celda
    pipeline.modifiers.append(AffineTransformationModifier(
        relative_mode=False,
        target_cell=pipeline.compute(0).cell[...]
    ))

    # Aplicar el modificador SelectTypeModifier para seleccionar todas las partículas del tipo 'Ox'
    pipeline.modifiers.append(SelectTypeModifier(
        operate_on="particles",
        property="Particle Type",
        types={'Ox', 'Cgr'}
    ))

    # Agregar el modificador AffineTransformationModifier para aplicar transformaciones absolutas a la celda nuevamente
    pipeline.modifiers.append(AffineTransformationModifier(
        relative_mode=False,
        target_cell=data.cell[...]
    ))

    # Agregar el modificador DeleteSelectedModifier para eliminar las partículas seleccionadas
    pipeline.modifiers.append(DeleteSelectedModifier())

    # Agregar el modificador ConstructSurfaceModifier para construir una superficie basada en el método AlphaShape
    pipeline.modifiers.append(ConstructSurfaceModifier(
        method=ConstructSurfaceModifier.Method.AlphaShape,
        radius=3.0,
        identify_regions=True
    ))

    data = pipeline.compute()


    empty_volume = data.attributes['ConstructSurfaceMesh.empty_volume']

    total_volume = (lx * ly * lz) * 1e-30  # Convertir de A^3 a m^3
    free_volume = empty_volume * 1e-30  # Convertir de A^3 a m^3

    return free_volume, total_volume
