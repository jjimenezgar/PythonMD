# PythonMD
Repositorio de Herramientas de Simulación Molecular en Python
Este repositorio contiene una colección de scripts en Python que se utilizan para realizar análisis y cálculos relacionados con la simulación molecular. Estas herramientas se pueden utilizar para analizar y comprender propiedades estructurales y dinámicas de sistemas moleculares.

Contenido del repositorio
El repositorio contiene los siguientes scripts:

1. Determinación del ECSA
El script ecsa.py calcula el Área Superficial Electroquímicamente Activa (ECSA) de una nanopartícula metálica. Utiliza la librería MDAnalysis de Python para analizar los átomos de la nanopartícula que están en contacto directo con el resto del sistema simulado. El resultado del cálculo proporciona la ECSA en m²/g y el porcentaje de utilización del catalizador en los diferentes sistemas.

2. Función de Distribución Radial (RDF)
El script rdf.py calcula la Función de Distribución Radial (RDF) para un sistema molecular dado. Utiliza la librería MDAnalysis de Python para calcular y guardar el RDF en un archivo de salida especificado. El RDF proporciona información sobre la distribución de partículas en el espacio, ofreciendo conocimientos sobre la estructura e interacciones del sistema.

3. Análisis estructural del agua
El script cluster_analysis.py utiliza la herramienta OVITO con su módulo de Python para realizar el análisis de clusters de agua en un sistema molecular. Identifica conjuntos de partículas conectadas y determina sus tamaños basados en un criterio de vecindad. Este análisis permite comprender la estructura y organización del agua en el sistema.

4. Determinación de estructuras cristalinas (CNA)
El script cna.py determina la estructura cristalina de nanopartículas utilizando el método de análisis de los primeros vecinos (CNA) de OVITO. Calcula las fracciones de las estructuras HCP, FCC, ICO y otras presentes en el sistema y genera un informe que incluye los resultados. Este análisis es útil para comprender la organización atómica de las nanopartículas.

5. Deformación atómica en nanopartículas
El script atomic_strain.py calcula la deformación atómica de un sistema utilizando la herramienta Atomic Strain de OVITO. Calcula el tensor de deformación y el gradiente de deformación a nivel atómico en cada partícula, a partir del movimiento relativo de sus vecinos. Proporciona información sobre las deformaciones locales en las nanopartículas.

6. Determinación del MSD y el coeficiente de difusión
El script MSD_Diffusion.py calcula el desplazamiento cuadrático medio (MSD) y grafica el coeficiente de difusión del agua en diferentes regiones a lo largo del eje Z en una simulación molecular. Utiliza la librería MDAnalysis para cargar la trayectoria de la simulación, calcular el MSD y obtener el coeficiente de difusión mediante una regresión lineal. También genera una visualización de la densidad del coeficiente de difusión.

7. Coordinación de átomos de oxígeno en simulaciones de dinámica molecular
Este repositorio contiene un script en Python llamado "cn_water.py" basado en la librería MDAnalysis. El script permite calcular y visualizar el número de coordinación de un conjunto de átomos de oxígeno en una simulación de dinámica molecular. Proporciona una herramienta útil para analizar la distribución de agua en esta región crítica.

8. Cálculo de densidad del agua en la región de TPB
Este repositorio también incluye otro script llamado "density_2D.py" que utiliza la librería MDAnalysis para calcular y visualizar la densidad media del agua en la región de TPB en sistemas de dinámica molecular. El script generará un mapa de calor que muestra la distribución de densidad. Esta herramienta es útil para analizar la distribución de agua en esta región crítica.

9. Solubilidad y permeaación de oxígeno e hidrógeno en la TPB
En este repositorio, también encontrará scripts relacionados con la solubilidad y permeación de oxígeno e hidrógeno en la TPB. Estos scripts permiten determinar las concentraciones de oxígeno e hidrógeno en sistemas con diferentes contenidos de ionómero, así como la concentración de gas presente en la superficie del metal.
