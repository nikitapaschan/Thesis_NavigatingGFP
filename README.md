# Master Thesis Repository

This repository contains all the materials necessary for reproducing the simulation studies and visualizations from my master's thesis to the title "Navigating the Garden of Forking Paths: Modelling and Addressing Researcher Degrees of Freedom".

## Folder Structure

- **code/**: Contains all R scripts used for simulations and visualizations.
- **datasets/**: Includes the main datasets generated during the simulation studies.
- **plots/**: Contains all the figures used in the thesis.



## Reproducing the Study

To reproduce the simulation studies and analyses, follow these steps:

1. **Install Required Packages**:
   - Run the script that installs and loads all necessary packages located in the `code/` folder.

2. **Run Main Functions**:
   - Execute the R script for the main simulation functions, which is located in the `code/` folder.

3. **Simulation of Max Z and Replication Studies**:
   - In the `Simulation_maxZ/` subfolder, you will find code to generate the initial simulation of maximum absolute Z-values (`maxZ`) and corresponding replication studies.
   - Each simulation setting is separated by the parameter $\mu$ (mean), allowing for a comparison between data with and without signal.

4. **Simulation of Researcher Degrees of Freedom (RDF)**:
   - The `Simulation_RDF/` subfolder contains the code for the simulation of researcher degrees of freedom, exploring different analysis strategies.

This structure facilitates easy navigation and reproduction of the simulation results used in this thesis.
