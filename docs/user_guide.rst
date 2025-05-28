User Guide
==========

This user guide provides detailed information on the various modules available in the CRISP package, along with performance tips to help optimize their usage.

CRISP is organized into two main subpackages:

- ``simulation_utility``: Tools for preparing and processing simulation data
- ``data_analysis``: Methods for analyzing simulation results

Simulation Utility Modules
===========================

1. Atomic Indices Classification
--------------------------------

**Module**: ``CRISP.simulation_utility.atomic_indices``

This module allows for the classification of atoms based on custom-defined indices, making it easier to analyze specific subsets within your simulation data.

**Key Features**:

- Customizable classification based on distance criteria
- Automatic generation of atom type indices
- Cutoff-based neighbor identification

**Example Usage**:

.. code-block:: python

    from CRISP.simulation_utility.atomic_indices import run_atom_indices
    
    cutoffs = {
        ("O", "H"): 1.2,
        ("Si", "O"): 1.8,
        ("Al", "Si"): 3.2,
        ("O", "O"): 3.0,
    }
    
    run_atom_indices("./wrapped_traj.traj", "./indices_new/", 
                     frame_index=2, cutoffs=cutoffs)

**Performance Tips**:

- Choose appropriate cutoff distances to balance accuracy and computational cost
- Use frame_index parameter to analyze representative frames instead of entire trajectories

2. Atomic Trajectory Visualization
----------------------------------

**Module**: ``CRISP.simulation_utility.atomic_traj_linemap``

Visualizes atomic trajectories in 3D to understand motion and structural changes.

**Key Features**:

- Interactive 3D visualization of atomic trajectories
- Customizable frame skipping and atom selection
- HTML output for web-based visualization

**Example Usage**:

.. code-block:: python

    from CRISP.simulation_utility.atomic_traj_linemap import plot_atomic_trajectory
    import numpy as np
    
    oxygen_indices = np.load("./indices_detailed/ex_fram_ox.npy")
    plot_atomic_trajectory(
        traj_path="./wrapped_traj.traj",
        selected_indices=oxygen_indices,
        output_path="oxygen_trajectories.html",
        frame_skip=10
    )

**Performance Tips**:

- Use frame_skip to reduce computational load for long trajectories
- Select specific atom indices rather than visualizing all atoms

3. Subsampling
--------------

**Module**: ``CRISP.simulation_utility.subsampling``

Extract representative structures from trajectories using Farthest Point Sampling with SOAP descriptors.

**Key Features**:

- SOAP-descriptor based diversity assessment
- Farthest Point Sampling for systematic subset selection
- Convergence visualization

**Example Usage**:

.. code-block:: python

    from CRISP.simulation_utility.subsampling import subsample
    
    all_frames = subsample(
        filename="./trajectory.traj",
        n_samples=30,
        index_type="all",
        file_format="traj",
        skip=10,
        plot_subsample=True
    )

**Performance Tips**:

- Use convergence plots to determine optimal subset size
- Consider computational budget vs. required accuracy when selecting n_samples

4. Error Analysis
-----------------

**Module**: ``CRISP.simulation_utility.error_analysis``

Perform statistical error analysis on time-correlated simulation data using autocorrelation function and block averaging methods.

**Key Features**:

- Autocorrelation function analysis with τ_int calculation
- Block averaging with convergence assessment
- Comprehensive statistical uncertainty estimation

**Example Usage**:

.. code-block:: python

    from CRISP.simulation_utility.error_analysis import error_analysis
    import numpy as np
    
    data = np.load("energy.npy")
    res = error_analysis(data)
    print(res["acf_results"])
    print(res["block_results"])

**Performance Tips**:

- Use autocorrelation analysis for statistically rigorous error estimation
- Compare ACF and block averaging results for validation

5. Interatomic Distance Calculation
-----------------------------------

**Module**: ``CRISP.simulation_utility.interatomic_distances``

Calculate and save distance matrices between atoms for further analysis.

**Key Features**:

- Periodic boundary condition handling
- Selective atom type analysis
- Distance matrix storage for reuse

**Example Usage**:

.. code-block:: python

    from CRISP.simulation_utility.interatomic_distances import distance_calculation, save_distance_matrices
    
    full_dms, sub_dms = distance_calculation("./wrapped_traj.traj", 
                                           frame_skip=10, index_type=["O"])
    save_distance_matrices(full_dms, sub_dms, ["O"], 
                          output_dir="distance_calculations")

**Performance Tips**:

- Use frame_skip to reduce computational load
- Focus on specific atom types rather than all atoms

Data Analysis Modules
=====================

1. Contact and Coordination Analysis
------------------------------------

**Module**: ``CRISP.data_analysis.contact_coordination``

Analyze coordination environments and dynamic contacts between atoms.

**Key Features**:

- Coordination number analysis with custom cutoffs
- Contact persistence analysis
- Interactive visualization with heatmaps

**Example Usage**:

.. code-block:: python

    from CRISP.data_analysis.contact_coordination import coordination, contacts
    
    # Coordination analysis
    cn = coordination("./wrapped_traj.traj", "O", ["O"], 
                     {('O', 'O'): 2.5}, skip_frames=10, plot_cn=True)
    
    # Contact analysis
    sub_dm, cal_contacts = contacts("./wrapped_traj.traj", "O", ["O"], 
                                   {('O', 'O'): 2.5}, plot_contacts=True)

**Performance Tips**:

- Use appropriate cutoffs based on chemical knowledge
- Skip frames for long trajectories to reduce computation time

2. Hydrogen Bond Analysis
-------------------------

**Module**: ``CRISP.data_analysis.h_bond``

Analyze hydrogen bond networks using geometric criteria based on angle and distance cutoffs.

**Key Features**:

- Geometric criteria (angle and distance based)
- Network visualization and correlation matrices
- Time series analysis of H-bond counts

**Example Usage**:

.. code-block:: python

    from CRISP.data_analysis.h_bond import hydrogen_bonds
    
    h_bonds_per_frame = hydrogen_bonds(
        traj_path="./wrapped_traj.traj",
        frame_skip=10,
        acceptor_atoms=["O"],
        angle_cutoff=120,
        h_bond_cutoff=2.5,
        bond_cutoff=1.6,
        plot_count=True,
        plot_heatmap=True,
        plot_graph_frame=True
    )

**Performance Tips**:

- Use relaxed cutoffs for initial screening, then refine parameters
- Enable specific plots based on analysis needs to save computation time

3. Radial Distribution Function (RDF)
-------------------------------------

**Module**: ``CRISP.data_analysis.prdf``

Calculate radial distribution functions for spatial relationship analysis.

**Key Features**:

- Partial RDF (PRDF) calculation for specific atom pairs
- Automatic binning and periodic boundary condition handling
- Interactive and static visualization options

**Example Usage**:

.. code-block:: python

    from CRISP.data_analysis.prdf import analyze_rdf
    
    result = analyze_rdf(
        use_prdf=True,
        rmax=6,
        traj_path="./trajectory.traj",
        nbins=50,
        frame_skip=100,
        atomic_indices=(pt_indices, pt_indices),
        create_plots=True
    )

**Performance Tips**:

- Adjust rmax to focus on relevant distance ranges
- Use frame_skip for long trajectories
- Limit analysis to essential atom pairs

4. Mean Square Displacement (MSD)
---------------------------------

**Module**: ``CRISP.data_analysis.msd``

Calculate and analyze Mean Square Displacement to study diffusion properties.

**Key Features**:

- Two-step analysis: MSD calculation and diffusion coefficient extraction
- Linear fitting with statistical error estimation
- Visual representation with fit lines

**Example Usage**:

.. code-block:: python

    from CRISP.data_analysis.msd import calculate_save_msd, analyze_from_csv
    
    # Step 1: Calculate MSD
    msd_values, msd_times = calculate_save_msd(
        traj_file="./trajectory.traj",
        timestep_value=5000.0,
        indices_file="./indices.npy",
        output_file="msd_results.csv"
    )
    
    # Step 2: Analyze diffusion
    D, error = analyze_from_csv("msd_results.csv", plot_msd=True)

**Performance Tips**:

- Use representative atom selection instead of all atoms
- Consider trajectory subsampling for very long simulations

5. Clustering Analysis
----------------------

**Module**: ``CRISP.data_analysis.clustering``

Identify atomic clusters using DBSCAN clustering algorithms.

**Key Features**:

- DBSCAN clustering with periodic boundary conditions
- Single-frame and trajectory-based analysis
- 3D visualization and cluster statistics

**Example Usage**:

.. code-block:: python

    from CRISP.data_analysis.clustering import analyze_frame, analyze_trajectory
    
    # Single frame analysis
    analyzer = analyze_frame(
        traj_path="./trajectory.traj",
        atom_indices="./indices.npy",
        threshold=2.6,
        min_samples=2,
        custom_frame_index=8100
    )
    result = analyzer.analyze_structure(output_dir="clustering_results")
    
    # Trajectory analysis
    analysis_results = analyze_trajectory(
        traj_path="./trajectory.traj",
        indices_path="./indices.npy",
        threshold=3.0,
        min_samples=2,
        frame_skip=100
    )

**Performance Tips**:

- Tune threshold and min_samples parameters based on physical understanding
- Use frame_skip for trajectory analysis to reduce computational load
- Consider using different thresholds to explore various connectivity definitions

6. Volumetric Density
---------------------

**Module**: ``CRISP.data_analysis.volumetric_atomic_density``

Create 3D volumetric density maps to visualize the spatial distribution of atoms throughout a trajectory.

**Key Features**:

- 3D density mapping with interactive visualizations
- Optional 2D projections (XY, YZ, XZ) for cross-sectional analysis
- Flexible thresholding (relative or absolute)
- Data export in NPZ format for further analysis
- Framework integration with selective atom visualization

**Example Usage**:

.. code-block:: python

    from CRISP.data_analysis.volumetric_atomic_density import create_density_map
    
    # Basic density map with projections
    fig = create_density_map(
        traj_path="./trajectory.traj",
        indices_path="./atom_indices.npy",
        frame_skip=100,
        threshold=0.05,
        opacity=0.8,
        show_projections=True,
        save_density=True,
        output_dir="./density_analysis",
        output_file="density_map.html"
    )
    
    # Comparative analysis with different projection modes
    projection_options = [True, False]
    for show_projections in projection_options:
        proj_text = "with_projections" if show_projections else "no_projections"
        create_density_map(
            traj_path="./trajectory.traj",
            indices_path="./indices.npy",
            frame_skip=100,
            threshold=0.0,
            opacity=0.8,
            show_projections=show_projections,
            save_projection_images=show_projections,
            output_dir="./density_comparison",
            output_file=f"density_{proj_text}.html"
        )

**Performance Tips**:

- Use appropriate frame_skip values to balance accuracy and computational efficiency
- Start with threshold=0.0 to see all density regions, then adjust for focus areas
- Use save_density=True to avoid recomputation when adjusting visualization parameters
- Consider omit_static_indices to focus on mobile species in complex systems
- Adjust nbins parameter based on system size and required resolution

**Applications**:

- Zeolite guest molecule distribution analysis
- Diffusion pathway identification in porous materials
- Adsorption site characterization
- Phase behavior studies in confined systems
- Structural flexibility analysis

General Performance Tips
========================

**Memory Management**:

- Use frame skipping for long trajectories
- Focus analysis on specific atom types or regions of interest
- Save intermediate results to avoid recomputation

**Parallel Processing**:

- Consider using multiple cores for trajectory analysis
- Batch process multiple systems or conditions
- Use efficient data structures (numpy arrays) for large datasets

**Parameter Optimization**:

- Start with reasonable default parameters based on chemical knowledge
- Use convergence plots and validation against experimental data
- Perform sensitivity analysis to understand parameter effects

Further Reading
===============

For more detailed information, including advanced usage and troubleshooting, please refer to the :doc:`developer_guide` and :doc:`tutorials`.
