CRISP
==============================
[//]: # (Badges)
[![GitHub Actions Build Status](https://github.com/REPLACE_WITH_OWNER_ACCOUNT/CRISP/workflows/CI/badge.svg)](https://github.com/REPLACE_WITH_OWNER_ACCOUNT/CRISP/actions?query=workflow%3ACI)
[![codecov](https://codecov.io/gh/REPLACE_WITH_OWNER_ACCOUNT/CRISP/branch/main/graph/badge.svg)](https://codecov.io/gh/REPLACE_WITH_OWNER_ACCOUNT/CRISP/branch/main)


<img src="https://github.com/Indranil17/CRISP_HOST/blob/main/crisp_logo2.png" width="300">

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

# Requirements
This package is built around the ASE (Atomic Simulation Environment) and thus requires the installation of ASE, available at: [ASE Installation Guide](https://wiki.fysik.dtu.dk/ase/install.html).

The clustering uses DBSCAN (Density-Based Spatial Clustering of Applications with Noise) implemented by scikit-learn, available at: [scikit-learn Installation](https://scikit-learn.org/stable/install.html).

For interactive 3D plots of clustering, CRISP utilizes the seaborn package. Ensure you have seaborn installed by following the instructions at: [Seaborn Installation](https://seaborn.pydata.org/installing.html).

# Installation
1. Now for the testing period, just pull the folder to the local drive, and then go inside the "CRISP" folder locally inside the terminal. 

2. In the same terminal type    
```
pip install .
```
, then you can run the package in your local system environment globally.

3. In a virtual environment as a standalone installation, please install the dependencies via the command line, type 
```
pip install joblib statsmodels IPython pandas scikit-learn plotly seaborn jupyter networkx
```

# Current Modules

To learn how to use each of the modules, please visit [examples](https://github.com/Indranil17/CRISP_HOST/tree/main/example).
Inside the example folder, check the [Notebook](https://github.com/Indranil17/CRISP_HOST/blob/main/example/CRISP_example.ipynb) to see them in use and also check the other folders of the example folder to see the outputs.
It has two sub-packages and ten modules. Please see the package's UML map. 
 
<img src="https://github.com/Indranil17/CRISP_HOST/blob/main/crisp_map.png" width="800">



# Acknowledgments
The package is from the (Nano)Materials modelling group, at Charles University.  
Email: lukas.grajciar@natur.cuni.cz  

<img src="https://github.com/Indranil17/CRISP_HOST/blob/main/group_logo.png" width="50">


### Contributors
If there is something you would like to see added to this package or if you would like to contribute,  

please email us at sahai@natur.cuni.cz or daniel.willimetz@natur.cuni.cz

