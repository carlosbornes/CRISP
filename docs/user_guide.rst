User Guide
==========

This user guide provides detailed information on the various modules
available in the CRISP package, along with performance tips to help optimize their usage.

Modules Overview
================

1. Atom Indices Classification
------------------------------

**Module**: CRISP.atom_indices

This module allows for the classification of atoms based on custom criteria,
such as cutoff distances between atom pairs. It is useful for distinguishing different regions 
or types of atoms within a structure.

**Key Features**:

- Customizable classification based on distance criteria.
- Ability to classify multiple sets of atoms within a single structure.

**Example Usage**:

.. code-block:: python

    from CRISP.atom_indices import classify_atoms
    atom_classes = classify_atoms(structure, cutoff=2.5)

**Performance Tips**:

- **Optimize Cutoff Selection**: Choosing an appropriate cutoff distance is crucial. \
  A smaller cutoff will reduce computation time but may miss important interactions, \
  while a larger cutoff increases accuracy but may be computationally expensive.
- **Parallel Processing**: If working with large structures, consider parallelizing \
  the classification process using Python’s multiprocessing library to speed up the computations.


2. Partial Radial Distribution Function (PRDF)
----------------------------------------------

**Module**: CRISP.prdf

This module calculates the Partial Radial Distribution Function (PRDF)
for a given structure, providing insights into spatial relationships between different atom types.

**Key Features**:

- Calculates PRDF for selected atom pairs.
- Option to compute dynamic PRDF over simulation trajectories.

**Example Usage**:

.. code-block:: python

    from CRISP.prdf import calculate_prdf
    prdf = calculate_prdf(structure, atom_pairs=[(1, 2)], r_max=10.0)

**Performance Tips**:

- **Limit PRDF Calculation to Essential Atom Pairs**: Only calculate the PRDF for atom \
  pairs of interest to save computational resources.
- **Adjust the r_max Parameter**: Reduce the r_max value if long-range interactions \
  are not essential to your analysis, which can significantly speed up the calculations.

3. Mean Square Displacement (MSD) Calculation
---------------------------------------------

**Module**: CRISP.msd_plot

The MSD module helps in analyzing the dynamics of atoms within a simulation by \
calculating their Mean Square Displacement over time.

**Key Features**:

- Supports customized atom selection for MSD calculations.
- Provides a visual plot of the MSD over time.

**Example Usage**:

.. code-block:: python

    from CRISP.msd_plot import calculate_msd
    msd = calculate_msd(trajectory, atom_indices=[1, 2, 3])

**Performance Tips**:

- **Select Representative Atoms**: Instead of calculating MSD for all atoms, \
  select a few representative atoms, especially in large systems, to reduce the computational load.
- **Trajectory Subsampling**: If the trajectory is very long, consider subsampling the \
  frames to perform MSD calculations on a reduced dataset without losing significant accuracy.

4. Hydrogen Bond Analysis
-------------------------

**Module**: CRISP.h_bond

This module calculates hydrogen bonds within a structure or trajectory, 
allowing you to understand the donor-acceptor interactions based on customizable cutoffs for angles and distances.

**Key Features**:

- Adjustable cutoffs for bond angles and distances.
- Analysis of hydrogen bonds over time in a trajectory.

**Example Usage**:

.. code-block:: python

    from CRISP.h_bond import analyze_h_bonds
    h_bonds = analyze_h_bonds(trajectory, angle_cutoff=30, distance_cutoff=2.5)

**Performance Tips**:

- **Use Relaxed Cutoffs for Initial Screening**: Start with more relaxed angle and \
  distance cutoffs to quickly screen for potential hydrogen bonds, then refine the \
  parameters in subsequent analyses.
- **Parallel Processing for Large Trajectories**: For long trajectories, distribute \
  the workload across multiple processors to improve performance.

5. Coordination Analysis
------------------------

**Module**: CRISP.coord

The Coordination module provides insights into the coordination environment of atoms\
in a structure by analyzing the number and type of neighbors around each atom.

**Key Features**:

- Customizable cutoff distances for neighbor detection.
- Suitable for analyzing local environments in complex materials.

**Example Usage**:

.. code-block:: python

    from CRISP.coord import calculate_coordination
    coord = calculate_coordination(structure, cutoff=3.0)

**Performance Tips**:

- **Focus on Key Atoms**: Limit the coordination analysis to atoms of particular interest\
  (e.g., active sites) rather than the entire structure to save time.
- **Use Efficient Data Structures**: Utilize numpy arrays or similar data structures to\
  handle large datasets efficiently during the analysis.

6. Clustering Frame Analysis
----------------------------

**Module**: CRISP.clustering_FrameAnalysis

This module applies clustering algorithms to identify groups of atoms within a \
single frame of a simulation, revealing patterns and structural features.

**Key Features**:

- Uses DBSCAN clustering with periodic boundary conditions.
- Visualizes clustering results with 3D plots.

**Example Usage**:

.. code-block:: python

    from CRISP.clustering_FrameAnalysis import perform_clustering
    clusters = perform_clustering(structure, eps=1.5, min_samples=5)

**Performance Tips**:

- **Tune Clustering Parameters**: Adjust eps and min_samples parameters carefully \
  to balance between detecting meaningful clusters and computational efficiency.
- **Reduce Dimensionality**: If applicable, perform dimensionality reduction \
  (e.g., PCA) before clustering to speed up the process.

7. Clustering Trajectory Analysis
---------------------------------

**Module**: CRISP.clustering_TrajAnalysis

Extending the clustering analysis to entire trajectories, this module helps identify \
dynamic structural changes and transitions over time.

**Key Features**:

- Tracks the evolution of clusters across frames.
- Provides tools for analyzing temporal patterns.

**Example Usage**:

.. code-block:: python

    from CRISP.clustering_TrajAnalysis import analyze_trajectory
    traj_clusters = analyze_trajectory(trajectory, eps=1.5, min_samples=5)

**Performance Tips**:

- **Frame Selection**: Perform clustering on a subset of frames first to identify important \
  regions of the trajectory before extending to the entire dataset.
- **Cluster Tracking**: Use efficient algorithms to track clusters across frames, especially in \
  long trajectories, to avoid memory and performance bottlenecks.

8. Atom Correlation Matrix
--------------------------

**Module**: CRISP.atom_correlation

This module calculates and visualizes the correlation matrix for atom pairs, 
based on their proximity across frames in a simulation.

**Key Features**:

- Computes correlation between atom pairs over time.
- Useful for identifying correlated movements or interactions.

**Example Usage**:

.. code-block:: python

    from CRISP.atom_correlation import calculate_correlation
    correlation_matrix = calculate_correlation(trajectory, cutoff=3.0)

**Performance Tips**:

- **Limit Atom Pair Selection**: Focus on atom pairs that are likely to exhibit \
  significant correlations to reduce computational complexity.
- **Efficient Memory Management**: Use sparse matrix representations if the \
  correlation matrix is large but sparsely populated to save memory and speed up computations.

Further Reading
===============

For more detailed information, including advanced usage and troubleshooting, please refer to the [Developer Guide](developer_guide.rst).
