Introductory Tutorials
======================

This section provides an overview of the basic functionalities offered by CRISP. \
Each tutorial introduces a specific feature of the package, guiding you through typical use cases.

First verify the installation by importing CRISP and its dependencies

.. code:: python

    import CRISP


Atom Indices
^^^^^^^^^^^^^^^^^^^
CRISP allows for the classification and visualization of atoms based on \
custom-defined indices, making it easier to analyze specific subsets within your simulation data.

1. **Classify Atom Indices**

   To classify atom indices, follow these steps:

   - Provide the path to the trajectory file (.traj format).
   - Specify the output folder where the classified atom indices will be saved.
   - Define the cutoffs to classify the atoms in a dictionary format.

   .. code:: python

       from CRISP.atom_indices import run_atom_indices

       # Path to your ASE trajectory file
       file_path = "./wrapped_traj.traj"

       # Output folder to save the indices
       output_folder = './indices_new/'

       # Define the cutoffs dictionary 
       cutoffs = {
           ("O", "H"): 1.2,
           ("Si", "O"): 1.8,
           ("Al", "Si"): 3.2,
           ("O", "O"): 3.0,
       }

       # Run the atom_indices function and save the results with a specific frame index
       run_atom_indices(file_path, output_folder, frame_index=2, cutoffs=cutoffs)

   **Output:**

   After running the above code, you'll receive output similar to the following, 
   which details the number of indices found for each atom type and where 
   the results are saved:

   .. code-block:: text

       Length of H indices: 120
       Length of Si indices: 168
       Length of Al indices: 24
       Length of O indices: 432
       Outputs saved.
       Saved cutoff indices for O-H to ./indices_new/cutoff/O-H_cutoff.csv
       Saved cutoff indices for Si-O to ./indices_new/cutoff/Si-O_cutoff.csv
       Saved cutoff indices for Al-Si to ./indices_new/cutoff/Al-Si_cutoff.csv
       Saved cutoff indices for O-O to ./indices_new/cutoff/O-O_cutoff.csv

2. **Visualize Atom Indices**

   Once you have classified the atom indices, you can visualize them to gain a better understanding of the spatial distribution and relationships within your molecular system.

   - Load the saved atom indices for further analysis or visualization.
   - Use the `visualize` function to create a visual representation of the atom indices for a specific frame.

   .. code:: python

       import numpy as np

       # Load the saved indices for Al atoms
       al_indices = np.load("./indices_new/Al_indices.npy")

       # Print the loaded indices
       print("Al_Indices:", al_indices)

       from CRISP.visualize_atom_indices import visualize

       # Visualize the atom indices for the specified frame
       visualize(file_path, frame_index=2)

   **Output:**

   The output from loading the indices will look like this:

   .. code-block:: text

       Al_Indices: [168 169 170 171 172 173 174 175 176 177 178 179 180 181 182 183 184 185
                    186 187 188 189 190 191]               


   The specified frame will also be visualized. 

.. image:: /images/introductory_tutorials/frame_visulation.png
   :alt: Frame Atoms 



Atom Coordination
^^^^^^^^^^^^^^^^^^^
CRISP provides tools to understand and explore the spatial
arrangement of atoms by analyzing coordination patterns within 
your simulation data. The `calculate_atom_coordination` function 
allows you to compute the coordination numbers of atoms based on 
custom-defined cutoffs.


1. **Calculate Atom Coordination**

   To calculate the coordination number of specific atoms in your trajectory:

   - Provide the path to the trajectory file (.traj format).
   - Specify the path to the saved indices of the atoms you are interested in.
   - Define the cutoffs for coordination analysis in a dictionary format.
   - Specify the atom type for which the coordination number will be calculated.

   .. code:: python

       from CRISP import cn_atom

       # Path to your ASE trajectory file
       file_path = "./wrapped_traj.traj"

       # Path to the saved indices of the atoms of interest
       indices_path = "./indices_new/O_indices.npy"

       # Output file to save the coordination data
       output_file = "coordination_data.pkl"

       # Define the cutoffs and atom coordination
       cutoffs = {("O", "O"): 3.0}
       atom_cn = "O"

       # Call the function with the file paths and parameters
       cn_atom.calculate_atom_coordination(file_path, indices_path, output_file, frame_skip=2, cutoffs=cutoffs, atom=atom_cn)

   **Output:**

   After running the code, you will obtain the coordination types for 
   the atoms of interest. Here is an example of the output:

   .. code-block:: text

       First entry of the coordination types list: 
       {192: 6, 193: 6, 194: 6, 195: 6, 196: 6, 197: 6, 198: 6, 199: 6, 
        200: 6, 201: 7, 202: 6, 203: 6, 204: 6, 205: 6, 206: 6, 207: 6, 
        208: 6, 209: 6, 210: 6, 211: 6, 212: 7, 213: 7, 214: 8, 215: 7, 
        216: 7, 217: 7, 218: 7, 219: 7, 220: 6, 221: 6, 222: 6, 223: 6, 
        224: 6, 225: 6, 226: 6, 227: 6, 228: 6, 229: 6, 230: 6, 231: 6, 
        232: 6, 233: 6, 234: 6, 235: 6, 236: 6, 237: 6, 238: 6, 239: 6, 
        240: 6, 241: 6}

   This output represents the coordination number of each atom of interest 
   within the specified cutoff distance. 
   The data is saved in a pickle file (`coordination_data.pkl`) 
   for further analysis or visualization.


2. **Visualize Atom Coordination**

   After calculating the coordination numbers, you can analyze and visualize the coordination distribution:

   - Use the path to the pickle file generated from the previous step.
   - Optionally specify a path to save the plot of the coordination distribution.

   .. code:: python

       from CRISP.cn_atom_analysis import main as analyze_and_plot

       # Define the file paths
       pickle_file = './coordination_data.pkl'
       plot_output_file = None  # or specify a path to save the plot, e.g., 'coordination_plot.png'

       # Call the analysis and plotting function
       analyze_and_plot(pickle_file, plot_output_file)

   **Output:**

   The analysis will provide an overview of the coordination distribution among the atoms, along with an image plot. Below is an example of the output:

   .. code-block:: text

       Overall Average percentage of 1-coordinated atoms: 0.53%
       Overall Average percentage of 2-coordinated atoms: 3.75%
       Overall Average percentage of 3-coordinated atoms: 5.24%
       Overall Average percentage of 4-coordinated atoms: 1.49%
       Overall Average percentage of 5-coordinated atoms: 1.60%
       Overall Average percentage of 6-coordinated atoms: 72.37%
       Overall Average percentage of 7-coordinated atoms: 14.65%
       Overall Average percentage of 8-coordinated atoms: 0.38%

   The results are visualized in a coordination plot, which can be saved as an image if a path is specified.

.. image:: /images/introductory_tutorials/atom_cn.png
   :alt: Coordination Number Plot

Hydrogen-Bonding
^^^^^^^^^^^^^^^^^^^

Hydrogen bond interactions between donor and acceptor atoms can be 
calculated and visualized in your simulation. Use the `analyze_hydrogen_bonds` 
function to perform the analysis. 

1. **Run Hydrogen Bond Analysis**

   To perform hydrogen-bond analysis, follow these steps:
   
   - A trajectory file (.traj format).
   - Indices files for donors, acceptors, and hydrogens (in `.npy` format).
   - Define parameters for your analysis such as distance cutoffs and angle cutoff.
   
   .. code:: python

      from CRISP.h_bond import analyze_hydrogen_bonds

      # Path to your ASE trajectory file
      traj_path = './wrapped_traj.traj'

      # Path to the indices of donor atoms
      donor_indices = './indices_detailed/ex_fram_ox.npy'

      # Path to the indices of acceptor atoms
      acceptor_indices = './indices_detailed/ex_fram_ox.npy'

      # Path to the indices of hydrogen atoms
      hydrogen_indices = './indices_detailed/wat_h.npy'

      # Output file to save the hydrogen bond data
      output_file = './hydrogen_bonds.csv'

      # Frame skip parameter
      frame_skip = 2

      # Distance cutoffs
      donor_acceptor_distance = 3.5
      donor_hydrogen_distance = 1.2

      # Angle cutoff in degrees
      angle_cutoff = 30.0

      # Call the analysis function
      analyze_hydrogen_bonds(traj_path, donor_indices, acceptor_indices, hydrogen_indices, output_file, frame_skip, donor_acceptor_distance, donor_hydrogen_distance, angle_cutoff)

   **Output:**

   After running the code, the hydrogen bond \
   information will be saved in a CSV file. This file will include details on the hydrogen bonds detected in the simulation frames according to the specified criteria.

   Example output:

   .. code-block:: text

      Frame, Donor Index, Acceptor Index, Hydrogen Index, Distance (Donor-Acceptor), Distance (Donor-Hydrogen), Angle (Donor-Hydrogen-Acceptor)
      1, 105, 210, 314, 3.45, 1.15, 28.4
      2, 106, 211, 315, 3.50, 1.20, 29.0
      ...

   This CSV file can be opened with standard data analysis \
   tools for further examination and visualization.

2. **Visualize Hydrogen Bond Data**

   After analyzing hydrogen bonds, you can create various plots to visualize 
   the hydrogen bond data. Use the following functions to generate and save plots.

   - **2D-weighted Hydrogen Bond Plot**

     The `h_bond_2d` function generates a 2D-weighted plot of the hydrogen bonds. 
     This plot helps visualize the spatial distribution of hydrogen bonds 
     in your simulation.

     .. code:: python

        from CRISP.h_bond_analysis import h_bond_2d

        # Path to the CSV file with hydrogen bond data
        hydrogen_bonds_csv = 'hydrogen_bonds.csv'

        # Generate 2D hydrogen bond plot
        h_bond_2d(hydrogen_bonds_csv)

     .. image:: /images/introductory_tutorials/hydrogen_bond_2d_plot.png
        :alt: 2D plot of hydrogen bonds

   - **Total Hydrogen Bonds Per Frame**

     The `plot_total_hydrogen_bonds_per_frame` function plots the total 
     number of hydrogen bonds detected per frame. 
     This helps in understanding the variation in hydrogen bonding 
     throughout the simulation.

     .. code:: python

        from CRISP.h_bond_analysis import plot_total_hydrogen_bonds_per_frame

        # Path to the CSV file with total hydrogen bonds per frame
        total_hydrogen_bonds_csv = 'total_hydrogen_bonds_per_frame.csv'

        # Plot total hydrogen bonds per frame
        plot_total_hydrogen_bonds_per_frame(total_hydrogen_bonds_csv)

     .. image:: /images/introductory_tutorials/total_hydrogen_bonds_per_frame.png
        :alt: Plot of total hydrogen bonds per frame

   - **Count Double Donor Hydrogen Bonds**

     The `count_double_donor_hydrogen_bonds` function counts the 
     double donor (DD) and single donor (SD) hydrogen bonds for each frame. 
     The results are saved to a text file.

     .. code:: python

        from CRISP.h_bond_analysis import count_double_donor_hydrogen_bonds, write_results_to_file

        # Path to the CSV file with hydrogen bond data
        hydrogen_bonds_csv = 'hydrogen_bonds.csv'
        
        # Path to save the results
        hydrogen_bond_results_txt = 'hydrogen_bond_results.txt'

        # Count DD and SD hydrogen bonds
        frame_counts = count_double_donor_hydrogen_bonds(hydrogen_bonds_csv)

        # Write the results to the text file
        write_results_to_file(hydrogen_bond_results_txt, frame_counts)

   - **Plot Hydrogen Bond Counts**

     The `plot_hydrogen_bond_counts` function creates a plot showing the counts of double donor hydrogen bonds for each frame.

     .. code:: python

        from CRISP.h_bond_analysis import plot_hydrogen_bond_counts

        # Path to the text file with hydrogen bond counts
        hydrogen_bond_results_txt = 'hydrogen_bond_results.txt'

        # Plot hydrogen bond counts
        plot_hydrogen_bond_counts(hydrogen_bond_results_txt)

     .. image:: /images/introductory_tutorials/hydrogen_bond_counts_plot.png
        :alt: Plot of hydrogen bond counts

3. **H-Bond Networks**

   For visualizing hydrogen-bond networks, 
   you can generate plots of hydrogen bond connectivity 
   for specific frames or over the entire trajectory. 
   This allows for detailed exploration of the hydrogen 
   bonding network structure in your system.

   - **Visualize H-Bond Networks:**

     .. code:: python

        from CRISP.h_bond_visualization import visualize_hydrogen_bonds

        # Example usage with average=True
        csv_file = './hydrogen_bonds.csv'
        wat_array_path = './indices_detailed/ex_fram_ox.npy'

        # For a single frame
        visualize_hydrogen_bonds(csv_file, wat_array_path, frame_index=1, average=False)

        # For a trajectory
        visualize_hydrogen_bonds(csv_file, wat_array_path, average=True)

   **Output:**

   The visualization produces two plots:

   - **Connectivity Matrix:**

    .. image:: ../images/introductory_tutorials/hb_connectivity_matrix.png
        :alt: Plot of h-bond connectivity matrix 
        :align: center


   - **NetworkX Plot:**

    .. image:: ../images/introductory_tutorials/hb_networkx_plot.png
        :alt: Plot of h-bond network 
        :align: center

   - **Node Size in the Graph:**

     - Represents the frequency of the indices during the simulation. 
     - A larger node size indicates a more persistent hydrogen bond.

   - **Edge Width between Nodes:**

     - Indicates the number of paired hydrogen bonds during the simulation.
     - A thicker edge width signifies a stable hydrogen bond between nodes/indices, while multiple thin edges suggest less stable hydrogen bonds with multiple indices.

   Overall, a bigger node size tells you that an index has a 
   more persistent hydrogen bond, and a thicker edge width indicates \
   stable hydrogen bonds with specific indices. Multiple thin edges may 
   reflect hydrogen bonds with various indices but less stability.


Radial Distribution Function (RDF)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Perform Radial Distribution Function (RDF) analysis to \
investigate the spatial relationships between atoms in your simulation.

1. **Run RDF Analysis**

   To perform RDF analysis, follow these steps:

   - Define the pairwise atom types you want to analyze.
   - Set parameters such as the maximum distance (`rmax`) and the number of bins (`nbins`).
   - Specify paths to the trajectory file and index files for the atoms of interest.
   - Optionally, provide a custom filename for the RDF output.

   .. code:: python

      from CRISP.prdf import analyze_rdf

      # Parameters for RDF analysis
      pairwise = ('O', 'Al')  # Pairwise atom symbols
      rmax = 12.0             # Maximum distance of RDF
      nbins = 100             # Number of bins to divide RDF
      use_prdf = True         # Flag to use PRDF for the indices
      traj_path = './wrapped_traj.traj'  # Path to the trajectory file
      wat_array_path = './indices_detailed/ex_fram_ox.npy'  # Path to the water indices array
      al_array_path = './indices_detailed/al_frw.npy'       # Path to the aluminum indices array
      output_filename = 'prdf_o_al_custom'  # Optional: specify a custom filename without extension

      # Call the analyze_rdf function with optional parameters
      x_data_all, y_data_all = analyze_rdf(pairwise, rmax, traj_path, wat_array_path, al_array_path, nbins, use_prdf, frame_skip=2, output_filename=output_filename) 

   **Output:**

   After running the RDF analysis, you will obtain the RDF data, 
   which can be visualized to understand the spatial distribution of atoms. 
   The results will be saved in a pickle file for plotting and animation.

2. **Visualize RDF**

   To visualize the RDF data, you can plot the RDF and create animations 
   to observe changes across simulation frames.

   - **Plot the RDF:**

     .. code:: python

        from CRISP.prdf_plot import plot_rdf_from_pickle

        # Plotting the RDF
        plot_rdf_from_pickle("./PRDF/prdf_o_al_custom.pkl")

     **Output:**

     The RDF plot will illustrate the radial distribution of atoms, 
     showing how the density of atoms varies as a function of distance.

    .. image:: ../images/introductory_tutorials/rdf_plot.png
        :alt: Plot PRDF
        :align: center

   - **Animate the RDF:**

     .. code:: python

        from CRISP.prdf_plot import animate_rdf_from_pickle

        # Create an animation of the RDF
        animate_rdf_from_pickle("./PRDF/prdf_o_al_custom.pkl")

     **Output:**

     The animation will display the changes in the RDF as frames progress, providing a dynamic view of the spatial relationships over time.

    .. image:: ../images/introductory_tutorials/rdf_animation.gif
        :alt: Animate PPRDF
        :align: center


Mean-Square Displacement (MSD)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Quantify atomic movements over time by calculating and 
visualizing the Mean-Square Displacement (MSD) of atoms. 
This analysis provides insights into the dynamics and diffusion 
characteristics of atoms in your simulation.

1. **Calculate MSD**

   To calculate the MSD, you will need to load the trajectory 
   and atom indices, and then use the `DiffusionCoefficient` class 
   from CRISP to compute the MSD and diffusion coefficient.

   - Provide the path to your trajectory file (.traj format).
   - Supply the path to the atom indices file (.npy format).
   - Set the timestep for the simulation. (important to keep in mind 
     the time stepsizes and no. of skips)
   
   .. code:: python

      from CRISP.msd_plot import DiffusionCoefficient
      from ase.io import read
      import numpy as np

      # Path to your ASE trajectory file
      traj_file = "./wrapped_traj.traj"

      # Path to the indices of the atoms
      atom_indices = np.load("./indices_detailed/ex_fram_ox.npy")

      # Time step in femtoseconds (0.5fs time stepsize, 1000 skips)
      timestep = 0.5 * 1000  

      # Read the trajectory
      traj = read(traj_file, ":")

      # Create an instance of DiffusionCoefficient and perform calculations
      diffusion_instance = DiffusionCoefficient(traj, timestep, atom_indices.tolist())
      diffusion_instance.calculate()
      diffusion_instance.plot()

   **Output:**

   After running the code, you will obtain the diffusion 
   coefficient and Standard Error (STD):

   Example output:

   .. code-block:: text

      Diffusion Coefficient (m^2/sec^-1): 4.935403012693005e-09
      Standard Error: 1.17e-04


2. **Visualize MSD**

   To visualize the MSD, use the provided plots to examine 
   the diffusion behavior of atoms in your simulation. The log-log plot helps to analyze the scaling relationship between MSD and time, while the diffusion calculation plot provides detailed insights into the diffusion coefficient.

   .. code:: python

      from CRISP.msd_plot import DiffusionCoefficient
      diffusion_instance.plot()

   **Outputs:**

   - **Log-Log Plot of MSD vs. Time**: 

     .. image:: /images/introductory_tutorials/msd_log_log_plot.png
        :alt: Log-Log Plot of MSD
        :align: center

   - **Diffusion Calculation Plot**:

     .. image:: /images/introductory_tutorials/diffusion_calculation_plot.png
        :alt: Diffusion Calculation Plot
        :align: center


Clustering
^^^^^^^^^^^^^^^^^^^

Use advanced clustering algorithms like DBSCAN to uncover 
patterns in atomic arrangements. This section demonstrates 
how to analyze and visualize clustering in your simulation data, 
including handling periodic boundary conditions.

1. **Perform Clustering Analysis**

   To analyze atomic arrangements and identify clusters, 
   follow these steps:

   - Provide the path to your trajectory file (.traj format).
   - Load the atom indices from a `.npy` file.
   - Set the clustering parameters such as distance threshold 
     and minimum samples.

   .. code:: python

      from CRISP.clustering_FrameAnalysis import StructureAnalyzer
      import numpy as np

      # Path to your ASE trajectory file
      traj_file = './wrapped_traj.traj'

      # Path to the atom indices file
      atom_indices = np.load('./indices_detailed/ex_fram_ox.npy')

      # Parameters for DBSCAN clustering
      threshold = 3.5  # Distance threshold for clustering
      min_samples = 2  # Minimum number of samples to form a cluster
      custom_frame_index = -1  # Analyze the last frame by default

      # Create an instance of StructureAnalyzer and perform the analysis
      analyzer = StructureAnalyzer(traj_file, atom_indices, threshold, min_samples, custom_frame_index=custom_frame_index)
      analyzer.analyze_structure()

   **Output:**

   After running the analysis, you will obtain clustering information, including:

   - **Number of Clusters (including noise)**
   - **Number of Noise Points (Cluster Indices with label -1)**
   - **Silhouette Score (including and excluding noise)**
   - **Cluster Indices and Average Cluster Size**

   Example output:

   .. code-block:: text

      Number of Clusters (including noise): 10
      Number of Noise (Cluster Indices with label -1): 3
      Silhouette Score (including noise): 0.32688493265485136
      Silhouette Score (excluding noise): 0.3973633101025257
      Cluster Indices: {0: array([576, 577, 579, 580, 585, 586, 589, 595, 601, 604, 606, 608, 612,
         613, 617, 623]), 1: array([578, 584, 607]), 2: array([581, 592, 620]), 3: array([583, 600, 616]), 4: array([587, 588, 593, 598, 603, 619, 621]), 5: array([590, 610, 614]), 6: array([591, 611, 615, 622]), 7: array([594, 597]), 8: array([599, 602, 605, 618])}
      Average Cluster Size: 5.0

2. **Visualize Clustering**

   To visualize the clustering results, use the interactive plotting capabilities provided by `StructureAnalyzer`. This will generate an interactive Plotly plot that displays the clusters and their characteristics.

   .. code:: python

      # Generate and display interactive clustering plot
      analyzer.analyze_structure()

   **Output:**

   The analysis will produce an interactive plot that visualizes 
   the clusters in your simulation. The plot allows you to 
   explore the clustering results interactively.

   Example plot:

   .. image:: /images/introductory_tutorials/clustering_interactive_plot.png
      :alt: Interactive Clustering Plot
      :align: center


Atom Correlation
^^^^^^^^^^^^^^^^^^^

Analyze and visualize the correlation between atom pairs to gain insights into their dynamic relationships over the course of the simulation. This section demonstrates how to calculate and visualize correlations between atoms.

1. **Perform Atom Correlation Analysis**

   To analyze the correlation between atom pairs, follow these 
   steps:

   - Provide the path to your trajectory file (.traj format).
   - Load the indices of the atoms to be analyzed.
   - Set parameters such as cutoff distance and frame skipping.

   .. code:: python

      from CRISP.atom_correlation import plot_correlation_matrix

      # Path to your ASE trajectory file
      file_path = './wrapped_traj.traj'

      # Paths to the indices of atom pairs
      atom1_indices_path = "./indices_detailed/ex_fram_ox.npy"
      atom2_indices_path = "./indices_detailed/ex_fram_ox.npy"

      # Call the function to compute and save correlation data
      plot_correlation_matrix(
          file_path=file_path, 
          atom1_indices_path=atom1_indices_path, 
          atom2_indices_path=atom2_indices_path, 
          cutoff=2.8, 
          frame_skip=1, 
          average=False, 
          output_dir='./outputs'
      )

   **Output:**

   After running the analysis, you will receive the following
   correlation information:

   - **Average number of positive correlations per frame**
   - **Total number of positive correlations**
   - **Average Correlation per indices pair**

   Example output:

   .. code-block:: text

      Average number of positive correlations per frame: 108.00
      Total number of positive correlations (from atom1_atom2_positive_correlations.csv): 108
      Average Correlation per indices pair: 1.00

2. **Visualize Atom Correlation**

   To visualize the correlation data, two plots are generated:

   - **Correlation Matrix**: Displays the correlation values between atom pairs.
   - **Clustermap**: Provides a hierarchical clustering view of the correlation matrix, highlighting patterns in the data.

   **Example Visualizations:**

   - **Correlation Matrix**:

     .. image:: /images/introductory_tutorials/correlation_matrix.png
        :alt: Correlation Matrix
        :align: center

   - **Clustermap**:

     .. image:: /images/introductory_tutorials/correlation_clustermap.png
        :alt: Correlation Clustermap
        :align: center


