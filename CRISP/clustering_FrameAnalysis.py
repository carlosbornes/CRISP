import argparse
import numpy as np
from ase.io import read, Trajectory
from scipy.spatial.distance import pdist, squareform
from sklearn.cluster import DBSCAN
from sklearn.metrics import silhouette_score
import plotly.graph_objects as go
from ase.visualize import view

class StructureAnalyzer:
    def __init__(self, traj_file, atom_indices, threshold, min_samples, custom_frame_index=None):
        self.traj_file = traj_file
        self.atom_indices = atom_indices
        self.threshold = threshold
        self.min_samples = min_samples
        self.custom_frame_index = custom_frame_index
        self.labels_with_periodic = None

    def read_custom_frame(self):
        traj = Trajectory(self.traj_file)
        if self.custom_frame_index is not None:
            custom_frame = traj[self.custom_frame_index]
        else:
            custom_frame = traj[-1]
        return custom_frame

    def find_clusters(self, X, metric='precomputed'):
        db = DBSCAN(eps=self.threshold, min_samples=self.min_samples, metric=metric).fit(X)
        return db.labels_

    def create_3d_scatter_plot(self, X, labels, title):
        fig_3d = go.Figure()
        for label in np.unique(labels):
            cluster_points = X[labels == label]
            fig_3d.add_trace(
                go.Scatter3d(
                    x=cluster_points[:, 0],
                    y=cluster_points[:, 1],
                    z=cluster_points[:, 2],
                    mode='markers',
                    marker=dict(size=3),
                    name=f'Cluster {label}'
                )
            )
        fig_3d.update_layout(
            title=title,
            scene=dict(xaxis=dict(title='X'), yaxis=dict(title='Y'), zaxis=dict(title='Z'))
        )
        return fig_3d

    def analyze_structure(self):
        frame = self.read_custom_frame()
        X = np.array([frame.positions[i] for i in self.atom_indices])
        cell = frame.get_cell()

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

        self.labels_with_periodic = self.find_clusters(square, metric='precomputed')

        silhouette_avg = silhouette_score(square, self.labels_with_periodic)

        non_noise_indices = np.where(self.labels_with_periodic != -1)
        filtered_data = square[non_noise_indices]

        silhouette_avg_filtered = silhouette_score(filtered_data, self.labels_with_periodic[non_noise_indices])

        fig_3d_with_periodic = self.create_3d_scatter_plot(
            X, self.labels_with_periodic, 'Interactive 3D Scatter Plot - With Periodic Boundaries'
        )

        fig_3d_with_periodic.show()

        # Create a mapping from cluster indices to original atom indices
        cluster_indices = {}
        cluster_sizes = {}
        cluster_to_original = {}

        for cluster_id in np.unique(self.labels_with_periodic):
            if cluster_id != -1:
                cluster_indices[cluster_id] = np.where(self.labels_with_periodic == cluster_id)[0]
                cluster_sizes[cluster_id] = len(cluster_indices[cluster_id])
                cluster_to_original[cluster_id] = self.atom_indices[cluster_indices[cluster_id]]

        avg_cluster_size = np.mean(list(cluster_sizes.values()))

        num_clusters_with_noise = len(np.unique(self.labels_with_periodic))   
        
        print(f"\nNumber of Clusters (including noise): {num_clusters_with_noise}")
        print(f"Number of Noise (Cluster Indices with label -1): {np.sum(self.labels_with_periodic == -1)}")
        print("Silhouette Score (including noise):", silhouette_avg)
        print("Silhouette Score (excluding noise):", silhouette_avg_filtered)
        print("Cluster Indices:", cluster_to_original)
        print("Average Cluster Size:", avg_cluster_size)
        
        return 
        
    def visualize_clusters(self, custom_indices, custom_atom):
        # Read the structure from the file
        structure = self.read_custom_frame()

        # Create a copy of the structure to avoid modifying the original structure
        structure_copy = structure.copy()

        # Change symbols for custom indices
        for idx in custom_indices:
            structure_copy.symbols[idx] = custom_atom

        # Visualize the modified structure
        view(structure_copy)

        
        
def main(args):
    trajectory_path = args.trajectory_path
    atom_indices_path = args.atom_indices_path
    threshold = args.threshold
    min_samples = args.min_samples
    custom_frame_index = args.custom_frame_index
    custom_indices = eval(args.custom_indices) if args.custom_indices else None
    custom_atom = args.custom_atom if args.custom_atom else "F"
    
    # Load atom indices
    atom_indices = np.load(atom_indices_path)

    # Analyze the specified frame or the last frame
    analyzer = StructureAnalyzer(trajectory_path, atom_indices, threshold, min_samples, custom_frame_index=custom_frame_index)
    analyzer.analyze_structure()
    
    if custom_indices:
        analyzer.visualize_clusters(custom_indices, custom_atom)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Analyze trajectory with DBSCAN clustering")
    parser.add_argument("trajectory_path", type=str, help="Path to trajectory file (ASE compatible format).")
    parser.add_argument("atom_indices_path", type=str, help="Path to numpy array file containing atom indices.")
    parser.add_argument("threshold", type=float, help="DBSCAN clustering threshold.")
    parser.add_argument("--min_samples", type=int, help="Minimum number of samples in a cluster for DBSCAN.")
    parser.add_argument('--custom_indices', type=str, help='Custom cluster indices for visualization (e.g., [20, 21, 22, 23, 24]).')
    parser.add_argument('--custom_atom', type=str, help='Atom type to replace the custom indices for visualization (default: "F").')
    parser.add_argument('--custom_frame_index', type=int, help='Specific frame number to analyze. If not provided, the last frame will be analyzed.')

    
    args = parser.parse_args()
    
    main(args)
