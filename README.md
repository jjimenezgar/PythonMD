# PythonMD :snake: :computer:	
Repository of Molecular Dynamics Simulations Tools in Python

This repository contains a collection of Python scripts used for analysis and calculations related to molecular simulation. These tools can be utilized to analyze and understand the structural and dynamic properties of molecular systems.

# Repository Contents
The repository includes the following scripts:

1. ECSA Determination
The script ecsa.py calculates the Electrochemically Active Surface Area (ECSA) of a metallic nanoparticle. It utilizes the MDAnalysis library in Python to analyze the atoms of the nanoparticle that are in direct contact with the rest of the simulated system. The calculation result provides the ECSA in m²/g and the catalyst utilization percentage in different systems.

<p align="center">
  <a href="https://github.com/jjimenezgar/PythonMD/blob/master/ecsa.py"><img src="https://github.com/jjimenezgar/PythonMD/blob/master/Images/active_area_a.png" alt="ECSA" width="200px" />
</p>

2. Radial Distribution Function (RDF)
The script rdf.py calculates the Radial Distribution Function (RDF) for a given molecular system. It uses the MDAnalysis library in Python to calculate and save the RDF in a specified output file. The RDF provides information about the particle distribution in space, offering insights into the system's structure and interactions.

<p align="center">
  <a href="https://github.com/jjimenezgar/PythonMD/blob/master/rdf.py"><img src="https://github.com/jjimenezgar/PythonMD/blob/master/Images/waterbulk_rdf%20(1).png" alt="RDF" width="200px" />
</p>

3. Water Structural Analysis
The script cluster_analysis.py utilizes the OVITO tool with its Python module to perform cluster analysis of water in a molecular system. It identifies sets of connected particles and determines their sizes based on a neighborhood criterion. This analysis helps in understanding the structure and organization of water in the system.

<p align="center">
  <a href="https://github.com/jjimenezgar/PythonMD/blob/master/cluster_analysis.py"><img src="https://github.com/jjimenezgar/PythonMD/blob/master/Images/water_cluster_a.png" alt="Cluster" width="200px" />
</p>

4. Common Neighbor Analysis (CNA)
The script cna.py determines the crystal structure of nanoparticles using the Common Neighbor Analysis (CNA) method in OVITO. It calculates the fractions of HCP, FCC, ICO, and other present structures in the system and generates a report that includes the results. This analysis is useful for understanding the atomic organization of nanoparticles.

<p align="center">
  <a href="https://github.com/jjimenezgar/PythonMD/blob/master/cna.py"><img src="https://github.com/jjimenezgar/PythonMD/blob/master/Images/cna_pt%20(1).png" alt="Cluster" width="200px"/>
</p>

5. Atomic Strain in Nanoparticles
The script atomic_strain.py calculates the atomic deformation of a system using the Atomic Strain tool in OVITO. It calculates the deformation tensor and deformation gradient at the atomic level in each particle, based on the relative motion of its neighbors. It provides information about local deformations in nanoparticles.

<p align="center">
  <a href="https://github.com/jjimenezgar/PythonMD/blob/master/atomic_strain.py"><img src="https://github.com/jjimenezgar/PythonMD/blob/master/Images/strain_a%20(1).png" alt="Cluster" width="200px" />
</p>

6. MSD and Diffusion Coefficient Determination
The script MSD_Diffusion.py calculates the Mean Squared Displacement (MSD) and plots the diffusion coefficient of water in different regions along the Z-axis in a molecular simulation. It uses the MDAnalysis library to load the simulation trajectory, calculate the MSD, and obtain the diffusion coefficient through linear regression. It also generates a visualization of the diffusion coefficient density.

<p align="center">
  <a href="https://github.com/jjimenezgar/PythonMD/blob/master/MSD_Diffusion.py"><img src="https://github.com/jjimenezgar/PythonMD/blob/master/Images/cv_msd.png" alt="Cluster" width="200px" />
</p>


7. Water Molecule Coordination in Molecular Dynamics Simulations
This repository contains a Python script named cn_water.py based on the MDAnalysis library. The script allows for the calculation and visualization of the coordination number of a set of oxygen atoms in a molecular dynamics simulation. It provides a useful tool for analyzing the water distribution in this critical region.

<p align="center">
  <a href="https://github.com/jjimenezgar/PythonMD/blob/master/cn_water.py"><img src="https://github.com/jjimenezgar/PythonMD/blob/master/Images/CN_apendice.png" alt="Cluster" width="200px" />
</p>


8. Water Density Calculation in the TPB Region
This repository also includes another script named density_2D.py that uses the MDAnalysis library to calculate and visualize the average water density in the TPB region in molecular dynamics systems. The script generates a heatmap that displays the density distribution. This tool is helpful for analyzing the water distribution in this critical region.

<p align="center">
  <a href="https://github.com/jjimenezgar/PythonMD/blob/master/density_2D.py"><img src="https://github.com/jjimenezgar/PythonMD/blob/master/Images/water_density.png" alt="Cluster" width="150px"/>
</p>

10. Oxygen and Hydrogen Solubility and Permeation in the TPB
In this repository, you will also find scripts related to the solubility and permeation of oxygen and hydrogen in the TPB. These scripts allow for determining the concentrations of oxygen and hydrogen in systems with different ionomer contents, as well as the gas concentration present on the metal surface.

<p align="center">
  <a href="https://github.com/jjimenezgar/PythonMD/blob/master/solubility.py"><img src="https://github.com/jjimenezgar/PythonMD/blob/master/Images/sol_per_a.png" alt="Cluster" width="150px" />
</p>

Please refer to the individual scripts for more detailed information on their usage and functionalities.
