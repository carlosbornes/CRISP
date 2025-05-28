RDF (Radial Distribution Function)
==================================

Radial Distribution Functions are crucial to represent the structural correlations in the AIMD trajectories accurately. We reconstructed the results of the original work by computing the RDFs for all three pair types (O–O, O–H, H–H) across the five meta-GGA functionals. Our CRISP results are presented alongside experimental references at 298K provided by the original work gathered via X-ray diffraction (for g_OO), interpolated at 298K and joint X-ray/neutron diffraction experiments (for g_OH and g_HH).

Basic RDF Analysis
------------------

.. code:: python

    from CRISP.data_analysis.prdf import analyze_rdf
    import numpy as np
    import os

    # Define parameters
    temperatures = [750, 1000, 1250]
    pt_cases = ["Pt3", "Pt5"]

    for pt in pt_cases:
        # Load platinum indices
        pt_indices = np.load(f'./CHA_Data_Analysis/{pt}/indices_new/Pt_indices.npy')
        atomic_indices = (pt_indices.tolist(), pt_indices.tolist())
        
        for temp in temperatures:
            # Set paths
            pt_lower = pt.lower()
            traj_file = f"./CHA/{pt}/CHA_{pt_lower}_t{temp}.traj"
            output_dir = f"./Data_supplementary_analysis/PRDF/{pt}/{temp}K"
            
            # Run RDF analysis
            result = analyze_rdf(
                use_prdf=True,
                rmax=6,
                traj_path=traj_file,
                nbins=50,
                frame_skip=100,
                output_filename=f"{pt_lower}_t{temp}_ptpt",
                atomic_indices=atomic_indices,
                output_dir=output_dir,
                create_plots=True
            )

Results and Comparison
----------------------

**O–O RDFs:** CRISP replicates the overall features reported by Villard et al., with clear functional differences. M06-L and M11-L understructure water (weaker first peak, higher first minimum), while MN12-L and MN15-L lead to stronger, slightly overstructured networks. revM06-L achieves intermediate structuring. CRISP-derived RDFs systematically shift O–O peaks slightly closer to experimental positions.

**O–H and H–H RDFs:** CRISP and reported values are almost similar for O-H and H-H peaks, with minor details such as revM06-L showing more pronounced second peak (around 4.0 Å) aligning with experimental trends.

**Visualization:**

.. image:: ../images/specific_tutorials/rdf/case_study_rdf_page.jpg
   :width: 800
   :alt: Radial distribution functions computed using CRISP

Summary
-------

CRISP demonstrates optimized binning protocols, automated handling of periodic boundary conditions and statistical averaging, providing accurate reproduction of experimental trends across different theoretical methods.
