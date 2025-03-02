CRISP
==============================
[//]: # (Badges)
[![GitHub Actions Build Status](https://github.com/REPLACE_WITH_OWNER_ACCOUNT/CRISP/workflows/CI/badge.svg)](https://github.com/REPLACE_WITH_OWNER_ACCOUNT/CRISP/actions?query=workflow%3ACI)
[![codecov](https://codecov.io/gh/REPLACE_WITH_OWNER_ACCOUNT/CRISP/branch/main/graph/badge.svg)](https://codecov.io/gh/REPLACE_WITH_OWNER_ACCOUNT/CRISP/branch/main)


<img src="https://github.com/Indranil17/CRISP/blob/main/example/DALL%C2%B7E_art.png" width="300">

# CRISP (Clustering, Radial, and other Insightful Simulation Post-processing)

# What's CRISP?
CRISP is a post-simulation analysis package built on the foundation of Atomic Simulation Environment (ASE). It is designed to unlock the secrets hidden within your molecular dynamics simulations, offering a comprehensive suite of tools for in-depth exploration and analysis.

### **CRISP's Features:**
1. **Clustering Mastery:**
   - Uncover patterns in atomic arrangements through advanced clustering algorithms (DBSCAN) with built-in periodicity awareness.
   - Identify and visualize distinct groups of atoms, shedding light on structural transitions and dynamic behaviours.

2. **Radial Proficiency:**
   - Harness the Partial Radial Distribution Function (PRDF) analysis to dissect spatial relationships between atoms.
   - Tailor your analysis by defining custom indices, also providing dynamic PRDF of a simulation, hence further insights into the distribution of atomic pairs.

3. **Insightful Post-Processing - Customizable MSD:**
   - Illuminate the dynamics of your simulations with the Mean Square Displacement (MSD) calculation, offering a quantitative measure of atom movements.
   - Customize your analysis by specifying custom atom indices, allowing you to focus on specific subsets of atoms and gain targeted insights into their dynamics.

4. **Hydrogen Bond Elegance:**
   - Dive into the world of molecular interactions with hydrogen bond calculations.
   - Understand the subtle but crucial interactions of Donor-acceptor with customised based on the varied cutoff for angles and distances.

5. **Coordination Analysis:**
   - Gain a deeper understanding of the spatial arrangement of atoms by exploring coordination patterns in your simulation.
   - With customised cutoff for atom pairs gain further insights into the material's structure and behaviour.

7. **Correlation of Atom-Pairs:**
   - Computes and visualizes the correlation matrix between pairs of atoms based on their indices.
   - It analyzes the trajectory data to determine how often pairs of atoms are within a specified cutoff distance across frames.
     
8. **Atom Indices classification:**
   - A customisable method to classify the atoms of the structure with possible cutoffs for multiple pairs.
    


# Requirements
This package is built around the ASE (Atomic Simulation Environment) and thus requires the installation of ASE, available at: [ASE Installation Guide](https://wiki.fysik.dtu.dk/ase/install.html).

The clustering uses DBSCAN (Density-Based Spatial Clustering of Applications with Noise) implemented by scikit-learn, available at: [scikit-learn Installation](https://scikit-learn.org/stable/install.html).

For interactive 3D plots of clustering, CRISP utilizes the seaborn package. Ensure you have seaborn installed by following the instructions at: [Seaborn Installation](https://seaborn.pydata.org/installing.html).

Other packages rquired are the following:
   - joblib 
   - statsmodels 
   - pandas 
   - plotly 
   - networkx

# Installation
To be updated

# Current Modules

To learn how to use each of the modules, please visit [examples](https://github.com/Indranil17/TEST/tree/main/example).
   - CRISP.atom_indices
   - CRISP.prdf
   - CRISP.msd_plot
   - CRISP.h_bond
   - CRISP.coord
   - CRISP.clustering_FrameAnalysis
   - CRISP.clustering_TrajAnalysis
   - CRISP.atom_correlation

# Acknowledgments
The package is from the (Nano)Materials modelling group, at Charles University.
Email: lukas.grajciar@natur.cuni.cz

### Contributors
If there is something you would like to see added to this package or if you would like to contribute, please email us at sahai@natur.cuni.cz or daniel.willimetz@natur.cuni.cz

### Testers

