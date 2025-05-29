RDF (Radial Distribution Function)
==================================

Radial Distribution Functions are crucial to represent the structural correlations in the AIMD trajectories accurately. We reconstructed the results of the original work by Villard et al. [1]_ by computing the RDFs for all three pair types (O–O, O–H, H–H) across the five meta-GGA functionals. Our CRISP results are presented alongside experimental references at 298K provided by the original work gathered via X-ray diffraction (for g_OO), interpolated at 298K and joint X-ray/neutron diffraction experiments (for g_OH and g_HH).

Basic RDF Analysis
------------------

.. code:: python

    from CRISP.data_analysis.prdf import analyze_rdf
    import numpy as np
    import os

    # Minnesota functionals from Villard et al. study
    functionals = ['M06-L', 'M11-L', 'MN12-L', 'MN15-L', 'revM06-L']

    # Directory mapping for file paths
    func_dirs = {
        'M06-L': 'M06L',
        'M11-L': 'M11L', 
        'MN12-L': 'MN12L',
        'MN15-L': 'MN15L',
        'revM06-L': 'REVM06L'
    }

    peak_data = {}

    print("Analyzing O-H radial distribution functions...")

    for functional in functionals:
        dir_name = func_dirs[functional]
        print(f"\nWorking on {functional}...")
        
        # Set up file paths
        oxygen_indices = f'./Data_supplementary_analysis/MetaGGA/{dir_name}/indices_new/O_indices.npy'
        hydrogen_indices = f'./Data_supplementary_analysis/MetaGGA/{dir_name}/indices_new/H_indices.npy'
        trajectory = f'./supplementary_data/MetaGGA/{dir_name}/TRAJEC.traj'
        results_dir = f'./Data_supplementary_analysis/MetaGGA/{dir_name}/goh'
        
        os.makedirs(results_dir, exist_ok=True)
        
        try:
            # Load atom indices
            o_atoms = np.load(oxygen_indices)
            h_atoms = np.load(hydrogen_indices)
            
            # Setup for O-H partial RDF calculation
            atom_pairs = (o_atoms.tolist(), h_atoms.tolist())
            
            # Calculate RDF
            rdf_results = analyze_rdf(
                use_prdf=True,
                rmax=6.0,
                traj_path=trajectory,
                nbins=50,
                frame_skip=1,
                output_filename="prdf_o-atoms_h-atoms",
                atomic_indices=atom_pairs,
                output_dir=results_dir,
                create_plots=True
            )
            
            # Extract first coordination shell peak
            distances = rdf_results['x_data']
            avg_rdf = np.mean(rdf_results['y_data_all'], axis=0)
            first_peak_idx = np.argmax(avg_rdf)
            first_peak_pos = distances[first_peak_idx]
            
            peak_data[functional] = first_peak_pos
            print(f"  First O-H peak: {first_peak_pos:.2f} Å")
            
        except Exception as error:
            print(f"  Failed to process {functional}: {error}")

    # Results summary
    print("\n" + "="*40)
    print("O-H COORDINATION ANALYSIS SUMMARY")
    print("="*40)
    print("Functional\t\tFirst Peak (Å)")
    print("-"*40)

    for func in functionals:
        if func in peak_data:
            tab = "\t\t" if len(func) <= 6 else "\t"
            print(f"{func}{tab}{peak_data[func]:.2f}")

    print(f"\nAnalyzed {len(peak_data)} functionals successfully")

Similarly, we calculated gOO and gHH by modifying the atom indices appropriately:

- **gOO analysis**: Used oxygen indices for both reference and target atoms
- **gHH analysis**: Used hydrogen indices for both reference and target atoms

Results and Comparison
----------------------

**O–O RDFs:** CRISP replicates the overall features reported by Villard et al. [1]_, with clear functional differences. M06-L and M11-L understructure water (weaker first peak, higher first minimum), while MN12-L and MN15-L lead to stronger, slightly overstructured networks. revM06-L achieves intermediate structuring. CRISP-derived RDFs systematically shift O–O peaks slightly closer to experimental positions.

**O–H and H–H RDFs:** CRISP and reported values are almost similar for O-H and H-H peaks, with minor details such as revM06-L showing more pronounced second peak (around 4.0 Å) aligning with experimental trends.

**Visualization:**

.. image:: ../images/specific_tutorials/rdf/case_study_rdf_page.jpg
   :width: 600
   :alt: Radial distribution functions computed using CRISP

Summary
-------

CRISP demonstrates optimized binning protocols, automated handling of periodic boundary conditions and statistical averaging, providing accurate reproduction of experimental trends across different theoretical methods.

References
----------

.. [1] Villard, Justin, Bircher, Martin P., and Rothlisberger, Ursula. "Structure and dynamics of liquid water from ab initio simulations: adding Minnesota density functionals to Jacob's ladder." *Chemical Science* 15, no. 12 (2024): 4434-4451.
