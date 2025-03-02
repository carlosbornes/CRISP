# CRISP/data_analysis/clustering.py

import argparse
import numpy as np
from ase.io import read
from sklearn.cluster import DBSCAN
from sklearn.metrics import silhouette_score
import plotly.graph_objects as go
import pickle
import csv
import matplotlib.pyplot as plt
import os

class StructureAnalyzer:
    def __init__(self, traj_file, atom_indices, threshold, min_samples, metric='precomputed', custom_frame_index=None):
        self.traj_file = traj_file
        self.atom_indices = atom_indices
        self.threshold = threshold
        self.min_samples = min_samples
        self.metric = metric
        self.custom_frame_index = custom_frame_index
        self.labels = None
        self.distance_matrix = None
        
    def read_custom_frame(self):
        """Read a specific frame or the last frame from the trajectory."""
        try:
            if self.custom_frame_index is not None:
                frame = read(self.traj_file, index=self.custom_frame_index)
            else:
                frame = read(self.traj_file, index='-1')
            return frame
        except Exception as e:
            print(f"Error reading trajectory: {e}")
            return None
            
    def calculate_distance_matrix(self, atoms):
        """Calculate a distance matrix with periodic boundary conditions."""
        positions = atoms.positions[self.atom_indices]
        
        if len(self.atom_indices) < self.min_samples:
            raise ValueError(f"Not enough atoms ({len(self.atom_indices)}) to form clusters with min_samples={self.min_samples}")
        
        full_dm = atoms.get_all_distances(mic=True)
        n_atoms = len(self.atom_indices)
        self.distance_matrix = np.zeros((n_atoms, n_atoms))
        
        for i, idx_i in enumerate(self.atom_indices):
            for j, idx_j in enumerate(self.atom_indices):
                self.distance_matrix[i, j] = full_dm[idx_i, idx_j]
        
        return self.distance_matrix, positions
        
    def find_clusters(self):
        """Find clusters using DBSCAN."""
        if self.distance_matrix is None:
            raise ValueError("Distance matrix must be calculated first")
            
        db = DBSCAN(
            eps=self.threshold, 
            min_samples=self.min_samples, 
            metric=self.metric
        ).fit(self.distance_matrix if self.metric == 'precomputed' else None)
        
        self.labels = db.labels_
        return self.labels

    def analyze_structure(self, save_html_path=None, output_dir=None):
        """Analyze the structure and find clusters.
        
        Parameters:
            save_html_path: Path to save HTML visualization (if None, auto-generated)
            output_dir: Directory to save all results (overrides path in save_html_path)
        """
        # Read the structure
        frame = self.read_custom_frame()
        if frame is None:
            return None
        
        # Set up output directory
        if output_dir is None and save_html_path is not None:
            output_dir = os.path.dirname(save_html_path)
        if not output_dir:
            output_dir = "clustering_results"
        
        # Create output directory if it doesn't exist
        os.makedirs(output_dir, exist_ok=True)
        print(f"\nSaving results to directory: {output_dir}")
        
        # Generate HTML save path if not provided
        if save_html_path is None:
            base_name = os.path.splitext(os.path.basename(self.traj_file))[0]
            save_html_path = os.path.join(output_dir, f"{base_name}_clusters.html")
            
        distance_matrix, positions = self.calculate_distance_matrix(frame)
        self.find_clusters()
        
        # Calculate silhouette score (excluding outliers)
        silhouette_avg = calculate_silhouette_score(self.distance_matrix, self.labels)

        # Create interactive 3D visualization
        create_html_visualization(
            positions=positions, 
            labels=self.labels,
            title='Interactive 3D Cluster Visualization', 
            save_path=save_html_path
        )

        # Extract cluster information
        cluster_info = extract_cluster_info(self.labels, self.atom_indices)
        num_clusters = cluster_info["num_clusters"]
        outlier_count = cluster_info["outlier_count"]
        avg_cluster_size = cluster_info["avg_cluster_size"]
        cluster_to_original = cluster_info["cluster_to_original"]

        # Print summary
        print_cluster_summary(num_clusters, outlier_count, silhouette_avg, avg_cluster_size, cluster_to_original)
        
        # Save detailed frame information to text file
        frame_info_path = os.path.join(output_dir, "frame_data.txt")
        save_frame_info_to_file(
            frame_info_path, 
            self.threshold, 
            self.min_samples,
            num_clusters, 
            outlier_count, 
            silhouette_avg, 
            avg_cluster_size, 
            cluster_to_original, 
            self.labels, 
            self.atom_indices
        )
        
        # Save analysis results as pickle file
        pickle_path = os.path.join(output_dir, "single_frame_analysis.pkl")
        result_data = {
            "num_clusters": num_clusters,
            "outlier_count": outlier_count,
            "silhouette_avg": silhouette_avg,
            "avg_cluster_size": avg_cluster_size,
            "cluster_to_original": cluster_to_original,
            "labels": self.labels,
            "positions": positions,
            "parameters": {
                "threshold": self.threshold,
                "min_samples": self.min_samples,
                "trajectory": self.traj_file,
                "frame_index": self.custom_frame_index
            }
        }
        
        with open(pickle_path, 'wb') as f:
            pickle.dump(result_data, f)
        print(f"Full analysis data saved to: {pickle_path}")
        
        return result_data

def create_html_visualization(positions, labels, title, save_path):
    """Create and save a 3D HTML visualization of clusters."""
    fig = go.Figure()
    
    # Add traces for each cluster
    for label in np.unique(labels):
        cluster_points = positions[labels == label]
        label_name = "Outliers" if label == -1 else f'Cluster {label}'
        marker_size = 5
        marker_color = 'gray' if label == -1 else None
        
        fig.add_trace(
            go.Scatter3d(
                x=cluster_points[:, 0],
                y=cluster_points[:, 1],
                z=cluster_points[:, 2],
                mode='markers',
                marker=dict(
                    size=marker_size,
                    color=marker_color
                ),
                name=label_name
            )
        )
    
    # Update layout
    fig.update_layout(
        title=title,
        scene=dict(
            xaxis=dict(title='X'),
            yaxis=dict(title='Y'),
            zaxis=dict(title='Z')
        ),
        legend=dict(itemsizing='constant')
    )
    
    # Save HTML file
    fig.write_html(save_path)
    print(f"3D visualization saved to {save_path}")

def calculate_silhouette_score(distance_matrix, labels):
    """Calculate silhouette score, handling edge cases."""
    try:
        non_outlier_mask = labels != -1
        if np.sum(non_outlier_mask) > 1:
            # Extract the sub-matrix for non-outlier points
            filtered_matrix = distance_matrix[np.ix_(non_outlier_mask, non_outlier_mask)]
            filtered_labels = labels[non_outlier_mask]
            return silhouette_score(filtered_matrix, filtered_labels, metric='precomputed')
        return 0
    except ValueError:
        return 0

def extract_cluster_info(labels, atom_indices):
    """Extract cluster information from labels."""
    cluster_indices = {}
    cluster_sizes = {}
    cluster_to_original = {}

    for cluster_id in np.unique(labels):
        if cluster_id != -1:  # Only count actual clusters (not outliers)
            cluster_indices[cluster_id] = np.where(labels == cluster_id)[0]
            cluster_sizes[cluster_id] = len(cluster_indices[cluster_id])
            cluster_to_original[cluster_id] = atom_indices[cluster_indices[cluster_id]]

    outlier_count = np.sum(labels == -1)
    num_clusters = len([label for label in np.unique(labels) if label != -1])
    
    # Calculate average cluster size
    avg_cluster_size = np.mean(list(cluster_sizes.values())) if cluster_sizes else 0

    return {
        "num_clusters": num_clusters,
        "outlier_count": outlier_count,
        "avg_cluster_size": avg_cluster_size,
        "cluster_sizes": cluster_sizes,
        "cluster_to_original": cluster_to_original
    }

def print_cluster_summary(num_clusters, outlier_count, silhouette_avg, avg_cluster_size, cluster_to_original):
    """Print a summary of clustering results."""
    print(f"\nNumber of Clusters: {num_clusters}")
    print(f"Number of Outliers: {outlier_count}")
    print(f"Silhouette Score: {silhouette_avg:.4f}")
    print(f"Average Cluster Size: {avg_cluster_size:.2f}")
    print("Cluster Information:")
    
    for cluster_id, atoms in cluster_to_original.items():
        print(f"  Cluster {cluster_id}: {len(atoms)} points")

def save_frame_info_to_file(file_path, threshold, min_samples, num_clusters, outlier_count, 
                           silhouette_avg, avg_cluster_size, cluster_to_original, labels, atom_indices):
    """Save detailed frame information to a text file."""
    with open(file_path, 'w') as f:
        f.write(f"DBSCAN Clustering Analysis Results\n")
        f.write(f"================================\n\n")
        f.write(f"Parameters:\n")
        f.write(f"  Threshold (eps): {threshold}\n")
        f.write(f"  Min Samples: {min_samples}\n\n")
        f.write(f"Results:\n")
        f.write(f"  Number of Clusters: {num_clusters}\n")
        f.write(f"  Number of Outliers: {outlier_count}\n")
        f.write(f"  Silhouette Score: {silhouette_avg:.4f}\n")
        f.write(f"  Average Cluster Size: {avg_cluster_size:.2f}\n\n")
        
        f.write(f"Detailed Cluster Information:\n")
        for cluster_id, indices in cluster_to_original.items():
            f.write(f"  Cluster {cluster_id}: {len(indices)} points\n")
            f.write(f"    Original atom indices: {indices.tolist()}\n\n")
        
        f.write(f"Outlier Information:\n")
        outlier_indices = atom_indices[labels == -1]
        f.write(f"  {len(outlier_indices)} outliers\n")
        if len(outlier_indices) > 0:
            f.write(f"    Original atom indices: {outlier_indices.tolist()}\n")
    
    print(f"Detailed frame data saved to: {file_path}")

def analyze_trajectory(trajectory_path, atom_indices_path, threshold, min_samples, skip_frames=10,
                      output_dir="contact_analysis", save_html_visualizations=True):
    """Analyze an entire trajectory with DBSCAN clustering."""
    atom_indices = np.load(atom_indices_path)
    os.makedirs(output_dir, exist_ok=True)
    
    try:
        trajectory = read(trajectory_path, index=f'::{skip_frames}')
        if not isinstance(trajectory, list):
            trajectory = [trajectory]
    except Exception as e:
        print(f"Error reading trajectory: {e}")
        return []
    
    num_frames = len(trajectory)
    results = []
    
    # Track the first and last frame data for visualization
    first_frame_data = None
    last_frame_data = None
    
    # Get base trajectory name for HTML filenames
    base_traj_name = os.path.splitext(os.path.basename(trajectory_path))[0]

    # Create a file to save per-frame information
    frame_data_file = os.path.join(output_dir, f"{base_traj_name}_frame_data.txt")
    with open(frame_data_file, 'w') as f_frames:
        f_frames.write(f"DBSCAN Clustering Analysis - Per-Frame Results\n")
        f_frames.write(f"==========================================\n\n")
        f_frames.write(f"Parameters:\n")
        f_frames.write(f"  Trajectory: {trajectory_path}\n")
        f_frames.write(f"  Threshold (eps): {threshold}\n")
        f_frames.write(f"  Min Samples: {min_samples}\n")
        f_frames.write(f"  Frame Skip: {skip_frames}\n\n")

        for frame_index, structure in enumerate(trajectory):
            frame_number = frame_index * skip_frames
            
            # Calculate distance matrix with periodic boundaries
            full_dm = structure.get_all_distances(mic=True)
            n_atoms = len(atom_indices)
            distance_matrix = np.zeros((n_atoms, n_atoms))
            for i, idx_i in enumerate(atom_indices):
                for j, idx_j in enumerate(atom_indices):
                    distance_matrix[i, j] = full_dm[idx_i, idx_j]
            
            # Extract positions (for visualization purposes)
            positions = structure.positions[atom_indices]
            
            # Find clusters
            db = DBSCAN(eps=threshold, min_samples=min_samples, metric='precomputed').fit(distance_matrix)
            labels = db.labels_
            
            # Store frame data for visualization
            frame_data = {
                'structure': structure.copy(),
                'positions': positions.copy(),
                'labels': labels.copy()
            }
            
            if frame_index == 0:
                first_frame_data = frame_data
            last_frame_data = frame_data
            
            # Calculate silhouette score
            silhouette_avg = calculate_silhouette_score(distance_matrix, labels)
            
            # Extract cluster information
            cluster_info = extract_cluster_info(labels, atom_indices)
            num_clusters = cluster_info["num_clusters"]
            outlier_count = cluster_info["outlier_count"]
            avg_cluster_size = cluster_info["avg_cluster_size"]
            cluster_indices = {k: np.where(labels == k)[0] for k in np.unique(labels) if k != -1}
            
            # Store results
            results.append((
                frame_number, 
                num_clusters,
                outlier_count,
                silhouette_avg, 
                avg_cluster_size
            ))
            
            # Write frame data to the file
            f_frames.write(f"Frame {frame_number}:\n")
            f_frames.write(f"  Number of Clusters: {num_clusters}\n")
            f_frames.write(f"  Number of Outliers: {outlier_count}\n")
            f_frames.write(f"  Silhouette Score: {silhouette_avg:.4f}\n")
            f_frames.write(f"  Average Cluster Size: {avg_cluster_size:.2f}\n")
            f_frames.write("  Detailed Cluster Information:\n")
            
            for label in sorted(cluster_indices.keys()):
                indices = cluster_indices[label]
                f_frames.write(f"    Cluster {label}: {len(indices)} points\n")
            
            f_frames.write("\n")
    
    print(f"Per-frame data saved to: {frame_data_file}")
    
    # Save HTML visualizations for first and last frames if requested
    if save_html_visualizations and first_frame_data and last_frame_data:
        # Create first frame visualization
        first_html_path = os.path.join(output_dir, f"{base_traj_name}_first_frame_clusters.html")
        create_html_visualization(
            positions=first_frame_data['positions'],
            labels=first_frame_data['labels'],
            title=f"Clusters in First Frame (Frame 0)",
            save_path=first_html_path
        )
        
        # Create last frame visualization
        last_frame_number = (len(results) - 1) * skip_frames
        last_html_path = os.path.join(output_dir, f"{base_traj_name}_last_frame_clusters.html")
        create_html_visualization(
            positions=last_frame_data['positions'],
            labels=last_frame_data['labels'],
            title=f"Clusters in Last Frame (Frame {last_frame_number})",
            save_path=last_html_path
        )
    
    return results

def save_analysis_results(analysis_results, output_dir="contact_analysis", output_prefix="clustering_results"):
    """Save analysis results to CSV, TXT, and PKL files in the specified output directory."""
    os.makedirs(output_dir, exist_ok=True)
    
    output_csv_file = os.path.join(output_dir, f"{output_prefix}.csv")
    output_txt_file = os.path.join(output_dir, f"{output_prefix}.txt")
    output_pickle_file = os.path.join(output_dir, f"{output_prefix}.pkl")

    # Save the analysis results to a CSV file
    with open(output_csv_file, 'w', newline='') as csvfile:
        csv_writer = csv.writer(csvfile)
        csv_writer.writerow([
            "Frame Number", 
            "Number of Clusters",
            "Number of Outliers",
            "Silhouette Score", 
            "Average Cluster Size"
        ])
        for result in analysis_results:
            csv_writer.writerow(result)

    # Save the analysis results to a text file
    with open(output_txt_file, 'w') as f:
        # Calculate averages
        frame_numbers = [result[0] for result in analysis_results]
        num_clusters = [result[1] for result in analysis_results]
        outlier_counts = [result[2] for result in analysis_results]
        silhouette_scores = [result[3] for result in analysis_results]
        avg_cluster_sizes = [result[4] for result in analysis_results]
        
        avg_num_clusters = np.mean(num_clusters)
        avg_outlier_count = np.mean(outlier_counts)
        avg_silhouette = np.mean(silhouette_scores)
        avg_cluster_size = np.mean(avg_cluster_sizes)
        
        f.write(f"DBSCAN Clustering Analysis Summary\n")
        f.write(f"================================\n\n")
        f.write(f"Average Values Across All Frames:\n")
        f.write(f"  Average Number of Clusters: {avg_num_clusters:.2f}\n")
        f.write(f"  Average Number of Outliers: {avg_outlier_count:.2f}\n")
        f.write(f"  Average Silhouette Score: {avg_silhouette:.4f}\n")
        f.write(f"  Average Cluster Size: {avg_cluster_size:.2f}\n\n")
        
        f.write(f"Analysis Results by Frame:\n")
        for result in analysis_results:
            frame_number, num_clusters, outlier_count, silhouette_avg, avg_cluster_size = result
            f.write(f"Frame {frame_number}:\n")
            f.write(f"  Number of Clusters: {num_clusters}\n")
            f.write(f"  Number of Outliers: {outlier_count}\n")
            f.write(f"  Silhouette Score: {silhouette_avg:.4f}\n")
            f.write(f"  Average Cluster Size: {avg_cluster_size:.2f}\n\n")

    # Save the analysis results to a pickle file
    with open(output_pickle_file, 'wb') as picklefile:
        pickle.dump(analysis_results, picklefile)

    print(f"Analysis results saved to directory: {output_dir}")
    
    return output_pickle_file

def plot_analysis_results(pickle_file, output_dir=None):
    """Plot analysis results from a pickle file and save to specified directory."""
    with open(pickle_file, 'rb') as f:
        analysis_results = pickle.load(f)

    # Extract data for plotting
    frame_numbers = [result[0] for result in analysis_results]
    num_clusters = [result[1] for result in analysis_results]
    outlier_counts = [result[2] for result in analysis_results]
    silhouette_scores = [result[3] for result in analysis_results]
    avg_cluster_sizes = [result[4] for result in analysis_results]
    
    # Calculate averages
    avg_num_clusters = np.mean(num_clusters)
    avg_outlier_count = np.mean(outlier_counts)
    avg_silhouette = np.mean(silhouette_scores)
    avg_avg_cluster_size = np.mean(avg_cluster_sizes)

    # Create subplots with individual sizes
    fig, axs = plt.subplots(4, 1, figsize=(18, 16), sharex=True)

    # Plot for avg_cluster_sizes
    axs[0].plot(frame_numbers, avg_cluster_sizes, color='red', linestyle='-', linewidth=2)
    axs[0].axhline(y=avg_avg_cluster_size, color='darkred', linestyle='--', alpha=0.7, 
                 label=f'Average: {avg_avg_cluster_size:.2f}')
    axs[0].set_ylabel('Average Cluster Size', color='red', fontsize=16)
    axs[0].tick_params(axis='y', labelcolor='red', labelsize=14)
    axs[0].grid(True, alpha=0.3)
    axs[0].legend(fontsize=12)

    # Plot for number of clusters
    axs[1].plot(frame_numbers, num_clusters, color='blue', linestyle='-', linewidth=2)
    axs[1].axhline(y=avg_num_clusters, color='darkblue', linestyle='--', alpha=0.7,
                 label=f'Average: {avg_num_clusters:.2f}')
    axs[1].set_ylabel('Number of Clusters', color='blue', fontsize=16)
    axs[1].tick_params(axis='y', labelcolor='blue', labelsize=14)
    axs[1].grid(True, alpha=0.3)
    axs[1].legend(fontsize=12)

    # Plot for outlier counts
    axs[2].plot(frame_numbers, outlier_counts, color='purple', linestyle='-', linewidth=2)
    axs[2].axhline(y=avg_outlier_count, color='darkviolet', linestyle='--', alpha=0.7,
                 label=f'Average: {avg_outlier_count:.2f}')
    axs[2].set_ylabel('Number of Outliers', color='purple', fontsize=16)
    axs[2].tick_params(axis='y', labelcolor='purple', labelsize=14)
    axs[2].grid(True, alpha=0.3)
    axs[2].legend(fontsize=12)

    # Plot for silhouette scores
    axs[3].plot(frame_numbers, silhouette_scores, color='orange', linestyle='-', linewidth=2)
    axs[3].axhline(y=avg_silhouette, color='darkorange', linestyle='--', alpha=0.7,
                 label=f'Average: {avg_silhouette:.4f}')
    axs[3].set_ylabel('Silhouette Score', color='orange', fontsize=16)
    axs[3].tick_params(axis='y', labelcolor='orange', labelsize=14)
    axs[3].grid(True, alpha=0.3)
    axs[3].legend(fontsize=12)
    axs[3].set_xlabel('Frame Number', fontsize=16)

    # Add a shared title
    plt.suptitle('Clustering Analysis Results', y=0.98, fontsize=20)

    # Adjust layout to prevent overlapping labels
    plt.tight_layout()
    plt.subplots_adjust(top=0.95)

    # Determine output directory if not provided
    if output_dir is None:
        output_dir = os.path.dirname(pickle_file)
    
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Create plot filename
    output_base = os.path.splitext(os.path.basename(pickle_file))[0]
    plot_file = os.path.join(output_dir, f"{output_base}_plot.png")
    
    # Save the figure
    plt.savefig(plot_file, dpi=300, bbox_inches='tight')
    
    # Show the plot
    plt.show()
    
    print(f"Analysis plot saved to: {plot_file}")

def main(args):
    trajectory_path = args.trajectory_path
    atom_indices_path = args.atom_indices_path
    threshold = args.threshold
    min_samples = args.min_samples
    output_dir = args.output_dir
    
    # Create base output directory
    os.makedirs(output_dir, exist_ok=True)
    print(f"Output files will be saved to: {output_dir}")
    
    # Load atom indices
    atom_indices = np.load(atom_indices_path)
    
    if args.mode == 'single':
        # Create a mode-specific subdirectory
        mode_dir = os.path.join(output_dir, "single_frame")
        os.makedirs(mode_dir, exist_ok=True)
        
        # Analyze a single frame
        analyzer = StructureAnalyzer(
            trajectory_path, 
            atom_indices, 
            threshold, 
            min_samples,
            metric='precomputed',
            custom_frame_index=args.custom_frame_index
        )
        
        # Define HTML output path
        base_traj_name = os.path.splitext(os.path.basename(trajectory_path))[0]
        frame_index = args.custom_frame_index if args.custom_frame_index is not None else "last"
        
        # Run analysis with explicit output directory
        analysis_result = analyzer.analyze_structure(output_dir=mode_dir)
        
    else:
        # Create a mode-specific subdirectory
        mode_dir = os.path.join(output_dir, "trajectory")
        os.makedirs(mode_dir, exist_ok=True)
        
        # Analyze the entire trajectory
        skip_frames = args.skip_frames
        
        # Perform trajectory analysis
        analysis_results = analyze_trajectory(
            trajectory_path, 
            atom_indices_path, 
            threshold, 
            min_samples, 
            skip_frames,
            output_dir=mode_dir,
            save_html_visualizations=True
        )
        
        # Save analysis results to the specified output directory
        pickle_file = save_analysis_results(
            analysis_results, 
            output_dir=mode_dir, 
            output_prefix=args.output_prefix
        )
        
        # Plot the analysis results
        plot_analysis_results(pickle_file, output_dir=mode_dir)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Analyze molecular structures with DBSCAN clustering."
    )
    parser.add_argument(
        "trajectory_path", 
        type=str, 
        help="Path to trajectory file (ASE compatible format)."
    )
    parser.add_argument(
        "atom_indices_path", 
        type=str, 
        help="Path to numpy array file containing atom indices."
    )
    parser.add_argument(
        "threshold", 
        type=float, 
        help="DBSCAN clustering threshold (eps parameter)."
    )
    parser.add_argument(
        "--min_samples", 
        type=int, 
        default=2, 
        help="Minimum number of samples in a cluster for DBSCAN. Default is 2."
    )
    parser.add_argument(
        "--mode", 
        type=str, 
        choices=["single", "trajectory"], 
        default="single", 
        help="Analysis mode: 'single' for single frame, 'trajectory' for whole trajectory."
    )
    parser.add_argument(
        "--output_dir", 
        type=str, 
        default="clustering", 
        help="Directory to save output files. Default is 'clustering'."
    )
    parser.add_argument(
        "--custom_frame_index", 
        type=int, 
        help="Specific frame number to analyze in 'single' mode. If not provided, the last frame will be analyzed."
    )
    parser.add_argument(
        "--skip_frames", 
        type=int, 
        default=10, 
        help="Skip frames in trajectory analysis. Default is 10."
    )
    parser.add_argument(
        "--output_prefix", 
        type=str, 
        default="clustering_results", 
        help="Prefix for output file names in trajectory analysis. Default is 'clustering_results'."
    )
    args = parser.parse_args()
    main(args)