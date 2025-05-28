from CRISP.data_analysis.clustering import analyze_trajectory, save_analysis_results, plot_analysis_results
import os
import numpy as np

# Define parameters for variations
temperatures = [750, 1000, 1250]
pt_cases = ["Pt3", "Pt5"]
threshold = 3.0
min_samples = 2
skip_frames = 100

for temp in temperatures:
    for pt in pt_cases:
        pt_lower = pt.lower()
        traj_file = f"./CHA/{pt}/CHA_{pt_lower}_t{temp}.traj"
        indices_file = f"./CHA_Data_Analysis/{pt}/indices_new/Pt_indices.npy"
        
        output_dir = f"{pt_lower}_{temp}K_traj_analysis_3.0"
        output_prefix = f"{pt_lower}_{temp}K_traj_clusters_3.0"
        
        print(f"Processing {pt} at {temp}K...")
        
        os.makedirs(output_dir, exist_ok=True)
        
        try:
            # Run trajectory analysis
            analysis_results = analyze_trajectory(
                traj_path=traj_file,
                indices_path=indices_file,
                threshold=threshold,
                min_samples=min_samples,
                frame_skip=skip_frames,
                output_dir=output_dir,
                save_html_visualizations=True  
            )
            
            pickle_file = save_analysis_results(
                analysis_results=analysis_results,
                output_dir=output_dir,
                output_prefix=output_prefix
            )
            
            plot_analysis_results(pickle_file, output_dir=output_dir)
            
