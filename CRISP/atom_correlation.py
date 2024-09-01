# correlation_matrix_plotter.py
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from ase.io import read
import csv
import pandas as pd
import os


def plot_correlation_matrix(
    file_path,
    atom1_indices_path,
    atom2_indices_path,
    cutoff=3.6,
    frame_skip=1,
    custom_frame=-1,
    average=False,
    output_dir="./",
):
    """
    Plot and analyze the correlation matrix between two types of atoms.

    Parameters:
        file_path (str): Path to the trajectory file.
        atom1_indices_path (str): Path to the file containing indices of the first type of atoms.
        atom2_indices_path (str): Path to the file containing indices of the second type of atoms.
        cutoff (float, optional): Cutoff distance to determine correlations. Default is 3.6.
        frame_skip (int, optional): Interval of frames to skip for averaging. Default is 1.
        custom_frame (int, optional): Specific frame number to analyze if not averaging. Default is -1 (last frame).
        average (bool, optional): Whether to average over all frames. Default is False.
        output_dir (str, optional): Directory to save the output files. Default is './'.

    Returns:
        None
    """
    # Load indices from provided paths
    atom1_indices = np.load(atom1_indices_path)
    atom2_indices = np.load(atom2_indices_path)

    num_atom1_indices = len(atom1_indices)
    num_atom2_indices = len(atom2_indices)
    correlation_sum = np.zeros((num_atom1_indices, num_atom2_indices))
    frame_count = 0
    positive_correlations_per_frame = []

    if average:
        traj = read(file_path, index=f"::{frame_skip}", format="traj")
        frame_count = len(traj)

        for frame in traj:
            dist_matrix = frame.get_all_distances(mic=True)
            selected_distances = dist_matrix[np.ix_(atom1_indices, atom2_indices)]
            correlation_matrix = (selected_distances < cutoff).astype(int)
            correlation_sum += correlation_matrix

            positive_correlations = []
            for i in range(num_atom1_indices):
                for j in range(num_atom2_indices):
                    if correlation_matrix[i, j] > 0:
                        positive_correlations.append(
                            (
                                (atom1_indices[i], atom2_indices[j]),
                                correlation_matrix[i, j],
                            )
                        )

            positive_correlations_per_frame.append(len(positive_correlations))

    else:
        frame = read(file_path, index=custom_frame, format="traj")
        dist_matrix = frame.get_all_distances(mic=True)
        selected_distances = dist_matrix[np.ix_(atom1_indices, atom2_indices)]
        correlation_matrix = (selected_distances < cutoff).astype(int)
        correlation_sum = correlation_matrix
        frame_count = 1

        positive_correlations = []
        for i in range(num_atom1_indices):
            for j in range(num_atom2_indices):
                if correlation_matrix[i, j] > 0:
                    positive_correlations.append(
                        ((atom1_indices[i], atom2_indices[j]), correlation_matrix[i, j])
                    )

        positive_correlations_per_frame.append(len(positive_correlations))

    average_correlation_matrix = correlation_sum / frame_count

    average_positive_correlations = sum(positive_correlations_per_frame) / frame_count
    print(
        f"Average number of positive correlations per frame: {average_positive_correlations:.2f}"
    )

    positive_correlations = []
    for i in range(num_atom1_indices):
        for j in range(num_atom2_indices):
            if average_correlation_matrix[i, j] > 0:
                positive_correlations.append(
                    (
                        (atom1_indices[i], atom2_indices[j]),
                        average_correlation_matrix[i, j],
                    )
                )

    os.makedirs(output_dir, exist_ok=True)

    csv_path_avg_pos_corr = os.path.join(
        output_dir, "atom1_atom2_average_positive_correlations_per_frame.csv"
    )
    with open(csv_path_avg_pos_corr, "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["Frame", "Average Positive Correlations"])
        for idx, count in enumerate(positive_correlations_per_frame):
            writer.writerow([idx + 1, count])

    txt_path_avg_pos_corr = os.path.join(
        output_dir, "atom1_atom2_average_positive_correlations_per_frame.txt"
    )
    with open(txt_path_avg_pos_corr, "w") as txtfile:
        for idx, count in enumerate(positive_correlations_per_frame):
            txtfile.write(f"Frame {idx + 1}: Average positive correlations = {count}\n")

    csv_path_pos_corr = os.path.join(
        output_dir, "atom1_atom2_positive_correlations.csv"
    )
    with open(csv_path_pos_corr, "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["Atom1 Index", "Atom2 Index", "Average Correlation"])
        for pair, value in positive_correlations:
            writer.writerow([pair[0], pair[1], f"{value:.2f}"])

    txt_path_pos_corr = os.path.join(
        output_dir, "atom1_atom2_positive_correlations.txt"
    )
    with open(txt_path_pos_corr, "w") as txtfile:
        for pair, value in positive_correlations:
            txtfile.write(
                f"Positive correlation between indices {pair[0]} and {pair[1]} with average correlation {value:.2f}\n"
            )

    print(
        f"Total number of positive correlations (from atom1_atom2_positive_correlations.csv): {len(positive_correlations)}"
    )

    # Read the CSV file into a DataFrame
    df = pd.read_csv(csv_path_pos_corr)

    # Calculate the average of the 'Average Correlation' column
    average_correlation = df["Average Correlation"].mean()

    # Print the average correlation
    print(f"Average Correlation per indices pair: {average_correlation:.2f}")

    fig, ax = plt.subplots(figsize=(10, 8))
    cax = ax.matshow(average_correlation_matrix, cmap="viridis")

    fig.colorbar(cax)

    ax.set_xticks(np.arange(num_atom2_indices))
    ax.set_yticks(np.arange(num_atom1_indices))

    ax.set_xticklabels(atom2_indices)
    ax.set_yticklabels(atom1_indices)

    plt.xticks(rotation=90)

    ax.set_xlabel("Atom2 Indices", fontsize=16)
    ax.set_ylabel("Atom1 Indices", fontsize=16)
    ax.set_title(
        "Average Correlation Matrix Between Atom1 and Atom2 Indices", fontsize=18
    )

    plt.tight_layout()
    plt.show()

    # Create a DataFrame from the correlation matrix for Seaborn
    df_corr = pd.DataFrame(
        average_correlation_matrix, index=atom1_indices, columns=atom2_indices
    )

    # Plot cluster map using Seaborn
    g = sns.clustermap(
        df_corr, method="average", cmap="RdBu", annot=True, annot_kws={"size": 8}
    )
    plt.setp(g.ax_heatmap.get_xticklabels(), rotation=60)
    plt.show()
