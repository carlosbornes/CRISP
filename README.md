CRISP
==============================
[//]: # (Badges)
[![GitHub Actions Build Status](https://github.com/REPLACE_WITH_OWNER_ACCOUNT/CRISP/workflows/CI/badge.svg)](https://github.com/REPLACE_WITH_OWNER_ACCOUNT/CRISP/actions?query=workflow%3ACI)
[![codecov](https://codecov.io/gh/REPLACE_WITH_OWNER_ACCOUNT/CRISP/branch/main/graph/badge.svg)](https://codecov.io/gh/REPLACE_WITH_OWNER_ACCOUNT/CRISP/branch/main)


<img src="https://github.com/Indranil17/CRISP/blob/main/example/DALL%C2%B7E_art.png" width="300">

# CRISP (Clustering, Radial, and other Insightful Simulation Post-processing)

## What is CRISP?
CRISP is a post-simulation analysis package built on the Atomic Simulation Environment (ASE). It is designed for efficient and insightful analysis of molecular dynamics (MD) and other simulations, enabling in-depth exploration with just a few lines of code, including powerful visualization options.

## Features
- **User-friendly**: Optimized for ease of use with detailed examples and extensive outputs for nuanced data analysis.
- **Highly parallelized**: Utilizes parallelization techniques that scale linearly with the number of CPU cores, allowing for fast analysis of large systems and long simulations on high-performance computing clusters.

## Analysis Toolkit

### 1. Cluster Comprehension
- Perform in-depth clustering analysis in molecular dynamics or Monte Carlo simulations using advanced algorithms like DBSCAN.
- Works with both periodic and non-periodic systems.
- Identify, visualize, and track distinct atom clusters to gain insights into unbiased clustering of selected atoms.

### 2. Customizable Radial Distribution Functions
- Compute and plot partial radial distribution functions (PRDF) for selected atoms or atom types.
- Easily analyze radial relationships between atoms with periodic boundary conditions.

### 3. Mean Square Displacement (MSD)
- Quantify atomic motion over time using MSD calculations, providing key insights into diffusion and dynamics.
- Customize analysis by selecting specific atom indices to focus on particular subsets of atoms.

### 4. Hydrogen Bond Analysis
- Identify and analyze hydrogen bonds with a single line of code.
- Customize hydrogen bond parameters or atom indices for detailed and specific analysis.
- Track structural parameters to understand the nature and stability of hydrogen bonds in your system.

### 5. Coordination Analysis
- Compute average coordination numbers for specific atom types with customizable cutoffs.
- Analyze contact times of selected atom types to study dynamic behavior efficiently.

### 6. Error Analysis
- Accurately estimate the error of any computed property using statistical techniques.
- Choose between autocorrelation function or block averaging to calculate the error of the mean, improving result reliability.
- Assess simulation convergence by analyzing vector or scalar properties like atomic positions or energy.

### 7. Efficient and Robust Sampling
- Sample structures using Furthest Point Sampling (FPS) with SOAP descriptors.
- Efficiently subsample large databases or simulations by selecting the most diverse structures while avoiding redundancy.

## Requirements
CRISP is primarily based on the Atomic Simulation Environment (ASE). You can install ASE following the guide: [ASE Installation Guide](https://wiki.fysik.dtu.dk/ase/install.html).

Additional dependencies include:
- `joblib`
- `statsmodels`
- `pandas`
- `plotly`
- `networkx`
- `DScribe`
- `seaborn`

## Installation
Installation instructions will be provided soon.
# Current Modules

To learn how to use each of the modules, please visit [examples]([https://github.com/Indranil17/TEST/tree/main/example](https://github.com/Indranil17/CRISP_HOST/blob/main/example/CRISP_latest_example.ipynb
)).
   - `CRISP.atom_indices`
   - `CRISP.prdf`
   - `CRISP.msd_plot`
   - `CRISP.h_bond`
   - `CRISP.coord`
   - `CRISP.clustering_FrameAnalysis`
   - `CRISP.clustering_TrajAnalysis`
   - `CRISP.atom_correlation`

# Acknowledgments
The package is from the (Nano)Materials modelling group, at Charles University.
Email: lukas.grajciar@natur.cuni.cz

### Contributors
If there is something you would like to see added to this package or if you would like to contribute, please email us at sahai@natur.cuni.cz or daniel.willimetz@natur.cuni.cz

### EXAMPLE
https://github.com/Indranil17/CRISP_HOST/blob/main/example/CRISP_latest_example.ipynb

### Docs
Master-Doc (Let's update in the Doc)
https://docs.google.com/document/d/14Bk8PkkE2FewRq1oMdXlqDkW0xUFwhyhYRJRxv559Ro/edit?usp=sharing

