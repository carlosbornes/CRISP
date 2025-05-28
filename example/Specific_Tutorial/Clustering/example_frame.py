from CRISP.data_analysis.clustering import analyze_frame
import os

# Paths and parameters
traj_file     = "./CHA/Pt5/CHA_pt5_t1250.traj"
indices_file  = "./CHA_Data_Analysis/Pt5/indices_new/Pt_indices.npy"
threshold     = 2.6
min_samples   = 2
frame_index   = 8100

# Create an output directory for this single‐frame analysis
output_dir = "clustering_results/Pt5_frames"
os.makedirs(output_dir, exist_ok=True)

# Initialize the analyzer
analyzer = analyze_frame(
    traj_path=traj_file,
    atom_indices=indices_file,
    threshold=threshold,
    min_samples=min_samples,
    custom_frame_index=frame_index
)

# Run the analysis and save all outputs (HTML, text, pickle)
result = analyzer.analyze_structure(
    save_html_path=os.path.join(output_dir, "frame8100_clusters.html"),
    output_dir=output_dir
)

# 'result' is a dict containing:
#   num_clusters, outlier_count, silhouette_avg, avg_cluster_size,
#   cluster_to_original, labels, positions, parameters
print("Single‐frame analysis complete. Results saved to:", output_dir)
print(result)