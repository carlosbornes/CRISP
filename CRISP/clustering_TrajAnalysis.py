import numpy as np
import argparse
from ase.io import Trajectory
from scipy.spatial.distance import pdist, squareform
from sklearn.cluster import DBSCAN
from sklearn.metrics import silhouette_score
import pickle
import csv
import matplotlib.pyplot as plt

def read_trajectory(file_path):
    return Trajectory(file_path, 'r')

def find_clusters(X, threshold, min_samples, metric='euclidean'):
    db = DBSCAN(eps=threshold, min_samples=min_samples, metric=metric).fit(X)
    return db.labels_

def calculate_jaccard_similarity(labels1, labels2):
    intersection = np.sum(labels1 == labels2)
    union = len(labels1) + len(labels2) - intersection
    return intersection / union

def analyze_trajectory(trajectory_path, atom_indices_path, threshold, min_samples, skip_frames=10):
    trajectory = read_trajectory(trajectory_path)
    atom_indices = np.load(atom_indices_path)
    
    num_frames = len(trajectory)
    results = []

    prev_labels = None  # Store cluster labels of the previous frame

    for frame_number, structure in enumerate(trajectory):
        if frame_number % skip_frames != 0:
            continue  # Skip this frame if it's not the 10th frame
        
        X = np.array([structure.positions[i] for i in atom_indices])
        cell = structure.get_cell()

        total = None
        for d in range(X.shape[1]):
            pd = pdist(X[:, d].reshape(X.shape[0], 1))
            pd[pd > np.linalg.norm(cell[d]) * 0.5] -= np.linalg.norm(cell[d])
            try:
                total += pd ** 2
            except Exception:
                total = pd ** 2

        total = np.sqrt(total)
        square = squareform(total)

        labels_with_periodic = find_clusters(square, threshold, min_samples, metric='precomputed')

        silhouette_avg = silhouette_score(X, labels_with_periodic)

        non_noise_indices = np.where(labels_with_periodic != -1)
        filtered_data = X[non_noise_indices]

        silhouette_avg_filtered = silhouette_score(filtered_data, labels_with_periodic[non_noise_indices])

        cluster_indices = {}
        cluster_sizes = {}

        for label in np.unique(labels_with_periodic):
            if label != -1:
                cluster_indices[label] = np.where(labels_with_periodic == label)[0]
                cluster_sizes[label] = len(cluster_indices[label])

        avg_cluster_size = np.mean(list(cluster_sizes.values()))

        # Calculate Jaccard similarity with previous frame's cluster assignments
        jaccard_similarity = 0.0
        if prev_labels is not None:
            jaccard_similarity = calculate_jaccard_similarity(prev_labels, labels_with_periodic)

        # Update previous labels for the next iteration
        prev_labels = labels_with_periodic.copy()

        results.append((frame_number, len(np.unique(labels_with_periodic)), silhouette_avg, silhouette_avg_filtered, cluster_indices, avg_cluster_size, jaccard_similarity))
    
    return results

def save_analysis_results(analysis_results, output_file_prefix):
    output_csv_file = f"{output_file_prefix}.csv"
    output_txt_file = f"{output_file_prefix}.txt"
    output_pickle_file = f"{output_file_prefix}.pkl"

    # Save the analysis results to a text (CSV) file
    with open(output_csv_file, 'w', newline='') as csvfile:
        csv_writer = csv.writer(csvfile)
        csv_writer.writerow(["Frame Number", "Number of Clusters", "Silhouette Score (including noise)", "Silhouette Score (excluding noise)", "Cluster Indices", "Average Cluster Size", "Jaccard Similarity"])
        for result in analysis_results:
            csv_writer.writerow(result)

    # Save the analysis results to a file
    with open(output_txt_file, 'w') as f:
        for result in analysis_results:
            frame_number, num_clusters, silhouette_avg_all, silhouette_avg_filtered, cluster_indices, avg_cluster_size, jaccard_similarity = result
            f.write(f"Frame {frame_number}:\n")
            f.write(f"Number of Clusters: {num_clusters}\n")
            f.write(f"Silhouette Score (including noise): {silhouette_avg_all}\n")
            f.write(f"Silhouette Score (excluding noise): {silhouette_avg_filtered}\n")
            
            # Cluster Indices
            f.write("Cluster Indices:\n")
            for label, indices in cluster_indices.items():
                f.write(f"  Cluster {label}: {indices}\n")
            
            f.write(f"Average Cluster Size: {avg_cluster_size}\n\n")
            f.write(f"Jaccard Similarity: {jaccard_similarity}\n\n")

    print(f"Analysis results saved to {output_txt_file}")

    # Save the analysis results to a pickle file
    with open(output_pickle_file, 'wb') as picklefile:
        pickle.dump(analysis_results, picklefile)

    print(f"Analysis results with changes saved to {output_csv_file} (CSV) and {output_pickle_file} (Pickle)")

def plot_analysis_results(pickle_file):
    # Load analysis results from the pickle file
    with open(pickle_file, 'rb') as f:
        analysis_results = pickle.load(f)

    # Extract data for plotting
    frame_numbers = [result[0] for result in analysis_results]
    num_clusters = [result[1] for result in analysis_results]
    silhouette_scores = [result[2] for result in analysis_results]
    silhouette_scores_exnoise = [result[3] for result in analysis_results]
    avg_cluster_sizes = [result[5] for result in analysis_results]
    jaccard_similarity = [result[6] for result in analysis_results]

    # Increase figure size before creating the subplots
    plt.figure(figsize=(20, 10))  # Adjust figure size for better visualization

    # Create four subplots with individual sizes
    fig, axs = plt.subplots(4, 1, figsize=(18, 16), sharex=True)

    # Plot for avg_cluster_sizes
    axs[0].plot(frame_numbers, avg_cluster_sizes, color='red', linestyle='-.')
    axs[0].set_ylabel('Avg Cluster Size', color='red', fontsize=18)
    axs[0].tick_params(axis='y', labelcolor='red', labelsize=15)
    axs[0].tick_params(axis='both', which='major', labelsize=15)

    # Set y-axis ticks and labels for better visualization
    yticks1 = axs[0].get_yticks()
    axs[0].set_yticks(yticks1[1:-1])
    axs[0].set_yticklabels([f'{int(label):,.0f}' for label in yticks1[1:-1]], fontsize=15)

    # Plot for num_clusters
    axs[1].plot(frame_numbers, num_clusters, color='green', linestyle=':')
    axs[1].set_ylabel('No. of Clusters', color='green', fontsize=18)
    axs[1].tick_params(axis='y', labelcolor='green', labelsize=15)
    axs[1].tick_params(axis='both', which='major', labelsize=15)

    # Set y-axis ticks and labels for better visualization
    yticks2 = axs[1].get_yticks()
    axs[1].set_yticks(yticks2[1:-1])
    axs[1].set_yticklabels([f'{int(label):,.0f}' for label in yticks2[1:-1]], fontsize=15)

    # Plot for silhouette_scores
    axs[2].plot(frame_numbers, silhouette_scores_exnoise, color='blue', linestyle='--')
    axs[2].set_ylabel('Silhouette Score', color='blue', fontsize=18)
    axs[2].tick_params(axis='y', labelcolor='blue', labelsize=15)
    axs[2].tick_params(axis='both', which='major', labelsize=15)

    # Set y-axis ticks and labels for better visualization
    yticks3 = axs[2].get_yticks()
    axs[2].set_yticks(yticks3[1:-1])
    axs[2].set_yticklabels([f'{label:.2f}' for label in yticks3[1:-1]], fontsize=15)

    # Plot for jaccard_similarity
    axs[3].plot(frame_numbers, jaccard_similarity, color='purple', linestyle=':')
    axs[3].set_xlabel('Frame Number', fontsize=18)
    axs[3].set_ylabel('Jaccard Similarity', color='purple', fontsize=16)
    axs[3].tick_params(axis='x', labelsize=15)
    axs[3].tick_params(axis='y', labelcolor='purple', labelsize=15)
    axs[3].tick_params(axis='both', which='major', labelsize=15)

    # Add a shared title
    plt.suptitle('Frame vs Analysis Results', y=1.02, fontsize=20)

    # Adjust layout to prevent overlapping labels
    plt.tight_layout()

    plt.show()

def main(args):
    trajectory_path = args.trajectory_path
    atom_indices_path = args.atom_indices_path
    threshold = args.threshold
    min_samples = args.min_samples
    skip_frames = args.skip_frames
    
    # Perform trajectory analysis
    analysis_results = analyze_trajectory(trajectory_path, atom_indices_path, threshold, min_samples, skip_frames)
    
    # Define output file prefix (common name for all output files)
    output_file_prefix = args.output_prefix
    
    # Save analysis results using the function
    save_analysis_results(analysis_results, output_file_prefix)
    
    # Plot the analysis results
    plot_analysis_results(f"{output_file_prefix}.pkl")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Analyze trajectory with DBSCAN clustering and calculate Jaccard similarity.")
    parser.add_argument("trajectory_path", type=str, help="Path to trajectory file (ASE compatible format).")
    parser.add_argument("atom_indices_path", type=str, help="Path to numpy array file containing atom indices.")
    parser.add_argument("threshold", type=float, help="DBSCAN clustering threshold.")
    parser.add_argument("--min_samples", type=int, help="Minimum number of samples in a cluster for DBSCAN.")
    parser.add_argument("--skip_frames", type=int, default=10, help="Skip frames in the analysis (default: 10).")
    parser.add_argument("--output_prefix", type=str, default="trajectory_analysis_results_with_changes", help="Prefix for output file names (default: trajectory_analysis_results_with_changes).")
    
    args = parser.parse_args()
    
    main(args)
