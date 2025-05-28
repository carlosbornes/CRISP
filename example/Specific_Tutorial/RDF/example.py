from CRISP.data_analysis.prdf import analyze_rdf
import numpy as np
import os

# Define parameters for variations
temperatures = [750, 1000, 1250]
pt_cases = ["Pt3", "Pt5"]

# Import the RDF analyzer
from CRISP.data_analysis.prdf import analyze_rdf

# Loop through all combinations
for pt in pt_cases:
    # Create main output directory for this Pt case
    main_output_dir = f"./Data_supplementary_analysis/PRDF/{pt}"
    os.makedirs(main_output_dir, exist_ok=True)
    
    # Load platinum indices
    pt_indices = np.load(f'./CHA_Data_Analysis/{pt}/indices_new/Pt_indices.npy')
    
    # For PRDF between platinum atoms, use the same indices for both reference and target
    atomic_indices = (pt_indices.tolist(), pt_indices.tolist())
    
    for temp in temperatures:
        print(f"\nProcessing {pt} at {temp}K...")
        
        # Construct trajectory path
        pt_lower = pt.lower()
        traj_file = f"./CHA/{pt}/CHA_{pt_lower}_t{temp}.traj"
        
        # Create temperature-specific output directory
        output_dir = f"{main_output_dir}/{temp}K"
        os.makedirs(output_dir, exist_ok=True)
        
        try:
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
            
            # Extract peak information
            x_data = result['x_data']
            y_data_avg = np.mean(result['y_data_all'], axis=0)
            peak_index = np.argmax(y_data_avg)
            peak_position = x_data[peak_index]
            
            print(f"  Pt-Pt coordination peak for {pt} at {temp}K: {peak_position:.2f} Å")
            print(f"  RDF data saved to: {output_dir}")
            