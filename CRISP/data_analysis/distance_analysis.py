"""
CRISP/data_analysis/distance_analysis.py

This module performs distance analysis on molecular dynamics trajectory data,
including coordination number calculation and atomic contact analysis.
"""

from ase.io import read
import numpy np
from joblib import Parallel, delayed
import pickle
from typing import Union, List, Dict, Tuple, Optional, Any
import os
import matplotlib.pyplot as plt
import itertools
from ase.data import vdw_radii, atomic_numbers, chemical_symbols
import seaborn as sns


def indices(atoms, ind: Union[str, List[Union[int, str]]]) -> np.ndarray:
    """
    Return array of atom indices from an ASE Atoms object based on the input specifier.
    
    Parameters
    ----------
    atoms : ase.Atoms
        ASE Atoms object containing atomic structure
    ind : Union[str, List[Union[int, str]]]
        Index specifier, can be "all", .npy file, integer(s), or chemical symbol(s)
        
    Returns
    -------
    np.ndarray
        Array of selected indices
        
    Raises
    ------
    ValueError
        If the index type is invalid
    """
    if ind == "all" or ind is None:
        return np.arange(len(atoms))
    if isinstance(ind, str) and ind.endswith(".npy"):
        return np.load(ind, allow_pickle=True)
    if not isinstance(ind, list):
        ind = [ind]
    if any(isinstance(item, int) for item in ind):
        return np.array(ind)
    if any(isinstance(item, str) for item in ind):
        idx = []
        if isinstance(ind, str):
            ind = [ind]
        for symbol in ind:
            idx.append(np.where(np.array(atoms.get_chemical_symbols()) == symbol)[0])
        return np.concatenate(idx)
    raise ValueError("Invalid index type")


def coordination_frame(atoms, target_atoms, ligand_atoms, custom_cutoffs=None, mic=True):
    """
    Calculate coordination numbers for target atoms based on interatomic distances and cutoff criteria.
    
    Parameters
    ----------
    atoms : ase.Atoms
        ASE Atoms object containing atomic structure
    target_atoms : Union[str, List[Union[int, str]]]
        Specifier for target atoms
    ligand_atoms : Union[str, List[Union[int, str]]]
        Specifier for ligand atoms
    custom_cutoffs : Optional[Dict[Tuple[str, str], float]]
        Dictionary with custom cutoff distances for atom pairs
    mic : bool
        Whether to use the minimum image convention
        
    Returns
    -------
    Dict[int, int]
        Dictionary mapping target atom indices to their coordination numbers
    """
    indices_target = indices(atoms, target_atoms)
    indices_ligand = indices(atoms, ligand_atoms)

    dm = atoms.get_all_distances(mic=mic)
    np.fill_diagonal(dm, np.inf)

    sub_dm = dm[indices_target, :][:, indices_ligand]

    target_atomic_numbers = np.array(atoms.get_atomic_numbers())[indices_target]
    ligand_atomic_numbers = np.array(atoms.get_atomic_numbers())[indices_ligand]

    target_vdw_radii = vdw_radii[target_atomic_numbers]
    ligand_vdw_radii = vdw_radii[ligand_atomic_numbers]

    cutoff_matrix = 0.6 * (target_vdw_radii[:, np.newaxis] + ligand_vdw_radii[np.newaxis, :])

    if custom_cutoffs is not None:
        cutoff_atomic_numbers = [tuple(sorted(atomic_numbers[symbol] for symbol in pair)) for pair in
                                 list(custom_cutoffs.keys())]
        cutoff_values = list(custom_cutoffs.values())

        cutoff_matrix_indices = [[tuple(sorted([i, j])) for j in ligand_atomic_numbers] for i in target_atomic_numbers]

        for i, target_atom in enumerate(cutoff_matrix_indices):
            for j, bond in enumerate(target_atom):
                if bond in cutoff_atomic_numbers:
                    cutoff_matrix[i, j] = cutoff_values[cutoff_atomic_numbers.index(bond)]

    coordination_numbers = np.sum(sub_dm < cutoff_matrix, axis=1)
    coordination_dict_frame = dict(zip(indices_target, coordination_numbers))
    return coordination_dict_frame


def get_avg_percentages(coordination_data, atom_type, plot=False, output_dir="./"):
    """
    Compute average percentages of each coordination number over all frames.
    
    Parameters
    ----------
    coordination_data : Dict[int, List[int]]
        Dictionary mapping atom indices to lists of coordination numbers
    atom_type : Optional[str]
        Chemical symbol of target atoms or None
    plot : bool, optional
        Boolean to indicate if a time series plot should be generated
    output_dir : str, optional
        Directory where output files will be saved
        
    Returns
    -------
    Dict[int, List[float]]
        Dictionary mapping each coordination number to a list of average percentages per frame
    """
    coord_types = set(itertools.chain.from_iterable(coordination_data.values()))
    coord_types = sorted(coord_types)

    num_frames = len(next(iter(coordination_data.values())))
    avg_percentages = {coord_type: [] for coord_type in coord_types}

    for frame_idx in range(num_frames):
        frame_data = [values[frame_idx] for values in coordination_data.values()]
        total_atoms = len(frame_data)
        for coord_type in coord_types:
            count = frame_data.count(coord_type)
            avg_percentage = count / total_atoms * 100
            avg_percentages[coord_type].append(avg_percentage)

    frames = list(range(len(next(iter(avg_percentages.values())))))
    coord_types = list(avg_percentages.keys())

    if plot:
        colors = plt.get_cmap('tab10', len(coord_types)).colors
        markers = itertools.cycle(['o', 's', 'D', '^', 'v', 'p', '*', '+', 'x'])
        plt.figure(figsize=(10, 6))
        for i, coord_type in enumerate(coord_types):
            plt.plot(frames, avg_percentages[coord_type], label=f'CN={coord_type}',
                     color=colors[i], marker=next(markers), markevery=max(1, len(frames) // 20))
        for i, coord_type in enumerate(coord_types):
            mean_value = sum(avg_percentages[coord_type]) / len(avg_percentages[coord_type])
            plt.axhline(y=mean_value, color=colors[i], linestyle='--', alpha=0.7,
                        label=f'Mean CN={coord_type}: {mean_value:.1f}%')
        plt.xlabel('Frame Index', fontsize=12)
        plt.ylabel('Percentage of Atoms (%)', fontsize=12)
        if atom_type is not None:
            plt.title(f'Coordination Analysis: {atom_type} Atoms', fontsize=14)
        else:
            plt.title('Coordination Analysis', fontsize=14)
        plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, "CN_time_series.png"), dpi=300, bbox_inches='tight')
        plt.close()

    return avg_percentages


def plot_coordination_distribution(avg_percentages, atom_type, plot, output_dir="./", output_file="CN_distribution"):
    """
    Plot a pie chart showing the overall average distribution of coordination numbers.
    
    Parameters
    ----------
    avg_percentages : Dict[int, List[float]]
        Dictionary of average percentages per coordination number
    atom_type : Optional[str]
        Chemical symbol for target atoms
    plot : bool
        Boolean to indicate if the plot should be generated
    output_dir : str, optional
        Directory where output files will be saved
    output_file : str, optional
        Filename for saving the plot
        
    Returns
    -------
    Dict[int, float]
        Dictionary of overall average coordination percentages
    """
    overall_avg_percentages = {coord_type: sum(percentages) / len(percentages) for coord_type, percentages in
                               avg_percentages.items()}

    if plot:
        plt.figure(figsize=(10, 7))
        labels = [f'CN={coord_type}: {pct:.1f}%' for coord_type, pct in overall_avg_percentages.items()]
        plt.pie(overall_avg_percentages.values(),
                labels=labels,
                autopct='%1.1f%%',
                colors=plt.cm.tab10.colors[:len(overall_avg_percentages)],
                startangle=90)
        plt.axis('equal')
        plt.title(f'Average Distribution of {atom_type} Atom Coordination', fontsize=14)
        plt.savefig(os.path.join(output_dir, f"{output_file}.png"), dpi=300, bbox_inches='tight')
        plt.close()

    return overall_avg_percentages


def log_coordination_data(distribution, atom_type, output_dir="./"):
    """
    Log coordination analysis statistics to a text file.
    
    Parameters
    ----------
    distribution : Dict[int, float]
        Dictionary of overall average coordination percentages
    atom_type : Optional[str]
        Chemical symbol of target atoms or None
    output_dir : str, optional
        Directory where the statistics file will be saved
        
    Returns
    -------
    None
    """
    if atom_type is not None:
        stats_file = f"CN_{atom_type}_statistics.txt"
    else:
        stats_file = "CN_statistics.txt"

    stats_file = os.path.join(output_dir, stats_file)
    with open(stats_file, 'w') as f:
        if atom_type is not None:
            f.write(f"Coordination Analysis for {atom_type} Atoms\n")
        else:
            f.write("Coordination Analysis\n")
        f.write("======================================\n\n")
        f.write("Overall Average Percentages:\n")
        for coord_type, avg_percentage in sorted(distribution.items()):
            f.write(f"  CN={coord_type}: {avg_percentage:.2f}%\n")


def coordination(filename, target_atoms, ligand_atoms, custom_cutoffs, skip_frames=10, 
                 plot=False, output_dir="./"):
    """
    Process a trajectory file to compute coordination numbers for each frame.
    
    Parameters
    ----------
    filename : str
        Path to the trajectory file
    target_atoms : Union[str, List[Union[int, str]]]
        Specifier for target atoms
    ligand_atoms : Union[str, List[Union[int, str]]]
        Specifier for ligand atoms
    custom_cutoffs : Dict[Tuple[str, str], float]
        Dictionary with custom cutoff distances
    skip_frames : int, optional
        Interval for skipping frames (default: 10)
    plot : bool, optional
        Boolean to indicate if plots should be generated (default: False)
    output_dir : str, optional
        Directory where output files will be saved (default: "./")
        
    Returns
    -------
    List
        List containing coordination dictionary, average percentages, distribution, and atom type
    """
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)

    trajectory = read(filename, index=f"::{skip_frames}")
    coordination_dict = {}

    # Process frames in parallel
    results = Parallel(n_jobs=-1)(
        delayed(coordination_frame)(atoms, target_atoms, ligand_atoms, custom_cutoffs)
        for atoms in trajectory
    )

    # Organize results by atom
    for frame_dict in results:
        for key, value in frame_dict.items():
            coordination_dict.setdefault(key, []).append(value)

    # Get atom type for labeling
    atom_type = (target_atoms if isinstance(target_atoms, str)
                 else chemical_symbols[target_atoms] if isinstance(target_atoms, int)
                 else None)

    # Generate analysis results and plots
    avg_percentages = get_avg_percentages(coordination_dict, atom_type, plot, output_dir=output_dir)
    distribution = plot_coordination_distribution(avg_percentages, atom_type, plot, output_dir=output_dir)
    log_coordination_data(distribution, atom_type, output_dir=output_dir)

    return [coordination_dict, avg_percentages, distribution, atom_type]


def contacts_frame(atoms, target_atoms, contact_atoms, custom_cutoffs, mic=True):
    """
    Processes a single atoms frame to compute the sub-distance matrix and the corresponding cutoff matrix.
    
    Parameters
    ----------
    atoms : ase.Atoms
        ASE Atoms object containing atomic structure
    target_atoms : Union[str, List[Union[int, str]]]
        Selection criteria for target atoms
    contact_atoms : Union[str, List[Union[int, str]]]
        Selection criteria for contact atoms
    custom_cutoffs : Dict[Tuple[str, str], float]
        Dictionary mapping atom pairs to custom cutoff values
    mic : bool, optional
        Whether to use minimum image convention (default: True)
        
    Returns
    -------
    Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]
        Sub-distance matrix, cutoff matrix, target indices, and contact indices
    """
    indices_target = indices(atoms, target_atoms)
    indices_contact = indices(atoms, contact_atoms)

    dm = atoms.get_all_distances(mic=mic)
    np.fill_diagonal(dm, np.inf)
    sub_dm = dm[np.ix_(indices_target, indices_contact)]

    target_atomic_numbers = np.array(atoms.get_atomic_numbers())[indices_target]
    contact_atomic_numbers = np.array(atoms.get_atomic_numbers())[indices_contact]

    target_vdw_radii = vdw_radii[target_atomic_numbers]
    contact_vdw_radii = vdw_radii[contact_atomic_numbers]

    cutoff_matrix = 0.6 * (target_vdw_radii[:, np.newaxis] + contact_vdw_radii[np.newaxis, :])

    if custom_cutoffs is not None:
        cutoff_atomic_numbers = [
            tuple(sorted(atomic_numbers[symbol] for symbol in pair))
            for pair in list(custom_cutoffs.keys())
        ]
        cutoff_values = list(custom_cutoffs.values())

        cutoff_matrix_indices = [
            [tuple(sorted([i, j])) for j in contact_atomic_numbers]
            for i in target_atomic_numbers
        ]
        for i, target_atom in enumerate(cutoff_matrix_indices):
            for j, bond in enumerate(target_atom):
                if bond in cutoff_atomic_numbers:
                    cutoff_matrix[i, j] = cutoff_values[cutoff_atomic_numbers.index(bond)]

    return sub_dm, cutoff_matrix, indices_target, indices_contact


def plot_contact_heatmap(contact_matrix, skip_frames, time_step, x_labels, y_labels, atom_type, output_dir="./"):
    """
    Plots and saves a heatmap showing contact times between target and contact atoms.
    
    Parameters
    ----------
    contact_matrix : np.ndarray
        Boolean 3D array of contacts
    skip_frames : int
        Number of frames skipped when processing the trajectory
    time_step : float
        Time step between frames (used to convert counts to time)
    x_labels : List[str]
        Labels for the contact atoms (x-axis)
    y_labels : List[str]
        Labels for the target atoms (y-axis)
    atom_type : Optional[str] 
        Target atom type (used for naming the file)
    output_dir : str, optional
        Directory where the output file will be saved (default: "./")
        
    Returns
    -------
    None
    """
    contact_times_matrix = np.sum(contact_matrix, axis=0) * skip_frames * time_step / 1000
    plt.figure(figsize=(10, 8))
    sns.heatmap(contact_times_matrix, xticklabels=x_labels, yticklabels=y_labels, 
                cmap="viridis", cbar_kws={"label": "Contact Time (ps)"})
    plt.xlabel("Contact Indices", fontsize=12)
    plt.ylabel("Target Indices", fontsize=12)
    plt.title("Heatmap of contact times within cutoffs", fontsize=14)
    plt.tight_layout()
    
    if atom_type is not None:
        filename = os.path.join(output_dir, f"{atom_type}_heatmap_contacts.png")
    else:
        filename = os.path.join(output_dir, "heatmap_contacts.png")
    
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    plt.close()


def plot_distance_heatmap(sub_dm_total, x_labels, y_labels, atom_type, output_dir="./"):
    """
    Plots and saves a heatmap of average distances between target and contact atoms.
    
    Parameters
    ----------
    sub_dm_total : np.ndarray
        3D numpy array containing sub-distance matrices for each frame
    x_labels : List[str]
        Labels for the contact atoms (x-axis)
    y_labels : List[str]
        Labels for the target atoms (y-axis)
    atom_type : Optional[str]
        Target atom type (used for naming the file)
    output_dir : str, optional
        Directory where the output file will be saved (default: "./")
        
    Returns
    -------
    None
    """
    average_distance_matrix = np.mean(sub_dm_total, axis=0)
    plt.figure(figsize=(10, 8))
    sns.heatmap(average_distance_matrix, xticklabels=x_labels, yticklabels=y_labels, 
                cmap="viridis", cbar_kws={"label": "Distance (Å)"})
    plt.xlabel("Contact Indices", fontsize=12)
    plt.ylabel("Target Indices", fontsize=12)
    plt.title("Distance matrix of selected atoms", fontsize=14)
    plt.tight_layout()
    
    if atom_type is not None:
        filename = os.path.join(output_dir, f"{atom_type}_heatmap_distance.png")
    else:
        filename = os.path.join(output_dir, "heatmap_distance.png")
    
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    plt.close()


def plot_contact_distance(sub_dm_total, contact_matrix, time_step, skip_frames, output_dir="./"):
    """
    Plots and saves a time series of the average contact distance over the trajectory.
    
    Parameters
    ----------
    sub_dm_total : np.ndarray
        3D numpy array of sub-distance matrices for each frame
    contact_matrix : np.ndarray
        Boolean 3D numpy array indicating which distances are considered contacts
    time_step : float
        Time between frames
    skip_frames : int
        Number of frames skipped when processing
    output_dir : str, optional
        Directory where the output file will be saved (default: "./")
        
    Returns
    -------
    None
    """
    contact_distance = np.where(contact_matrix, sub_dm_total, np.nan)
    contact_count = np.sum(contact_matrix, axis=(1, 2))/np.shape(sub_dm_total)[1]
    average_distance_contacts = np.nanmean(contact_distance, axis=(1, 2))

    x = np.arange(len(average_distance_contacts)) * time_step * skip_frames / 1000
    valid_indices = ~np.isnan(average_distance_contacts)
    interpolated = np.interp(x, x[valid_indices], average_distance_contacts[valid_indices])

    fig, ax1 = plt.subplots(figsize=(8, 6))

    ax1.plot(x, interpolated, '-', color="blue", label="Average Distance")
    ax1.scatter(x[valid_indices], average_distance_contacts[valid_indices], color="blue")
    ax1.axhline(np.mean(interpolated), linestyle="--", color="blue", 
                label=f"Mean: {np.mean(interpolated):.2f}")
    ax1.set_xlabel("Time [ps]", fontsize=12)
    ax1.set_ylabel("Distance [Å]", fontsize=12, color="blue")
    ax1.tick_params(axis='y', labelcolor="blue")

    ax2 = ax1.twinx()
    ax2.plot(x, contact_count, '-', color="red", label="Average Contact Count")
    ax2.axhline(np.mean(contact_count), linestyle="--", color="red", 
                label=f"Mean: {np.mean(contact_count):.1f}")
    ax2.set_ylabel("Contact count per atom", fontsize=12, color="red")
    ax2.tick_params(axis='y', labelcolor="red")

    plt.title("Average Distance of Contacts & Contact Count", fontsize=14)
    fig.tight_layout()

    fig.legend(loc="center right", bbox_to_anchor=(1.30, 0.5))

    filename = os.path.join(output_dir, "average_contact_analysis.png")
    plt.savefig(filename, dpi=300, bbox_inches="tight")
    plt.close()


def save_matrix_data(sub_dm_total, contact_matrix, output_dir="./"):
    """
    Saves the sub-distance matrices and contact matrix to npy files.
    
    Parameters
    ----------
    sub_dm_total : np.ndarray
        3D numpy array of sub-distance matrices
    contact_matrix : np.ndarray
        Boolean 3D numpy array of contact information
    output_dir : str, optional
        Directory where the output files will be saved (default: "./")
        
    Returns
    -------
    None
    """
    np.save(os.path.join(output_dir, "sub_dm_total.npy"), sub_dm_total)
    np.save(os.path.join(output_dir, "contact_matrix.npy"), contact_matrix)


def contacts(filename, target_atoms, contact_atoms, custom_cutoffs, skip_frames=10,
             plot_distance_matrix=False, plot_contacts=False, time_step=None, save_data=False,
             output_dir="./", mic=True):
    """
    Processes a molecular trajectory file to compute contacts between target and contact atoms.
    
    Parameters
    ----------
    filename : str
        Path to the trajectory file
    target_atoms : Union[str, List[Union[int, str]]]
        Criteria for selecting target atoms
    contact_atoms : Union[str, List[Union[int, str]]]
        Criteria for selecting contact atoms
    custom_cutoffs : Dict[Tuple[str, str], float]
        Dictionary with custom cutoff values for specific atom pairs
    skip_frames : int, optional
        Number of frames to skip (default: 10)
    plot_distance_matrix : bool, optional
        Boolean flag to plot average distance heatmap (default: False)
    plot_contacts : bool, optional
        Boolean flag to plot contact times heatmap (default: False)
    time_step : float, optional
        Time between frames in fs (required for contact heatmap)
    save_data : bool, optional
        Boolean flag to save matrices as npy files (default: False)
    output_dir : str, optional
        Directory where output files will be saved (default: "./")
    mic : bool, optional
        Whether to use minimum image convention (default: True)
        
    Returns
    -------
    Tuple[np.ndarray, np.ndarray]
        3D numpy array of sub-distance matrices and Boolean 3D numpy array of contacts
    """
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)

    sub_dm_list = []
    trajectory = read(filename, index=f"::{skip_frames}")

    # Process frames in parallel
    results = Parallel(n_jobs=-1)(
        delayed(contacts_frame)(atoms, target_atoms, contact_atoms, custom_cutoffs, mic=mic)
        for atoms in trajectory
    )

    sub_dm_list, cutoff_matrices, indices_target, indices_contact = zip(*results)
    cutoff_matrix = cutoff_matrices[0]
    sub_dm_total = np.array(sub_dm_list)

    # Get atom type for labeling
    atom_type = (target_atoms if isinstance(target_atoms, str)
                 else chemical_symbols[target_atoms] if isinstance(target_atoms, int)
                 else None)
                 
    y_labels = [f"{trajectory[0].get_chemical_symbols()[i]}({i})" for i in indices_target[0]]
    x_labels = [f"{trajectory[0].get_chemical_symbols()[i]}({i})" for i in indices_contact[0]]

    contact_matrix = sub_dm_total < cutoff_matrix

    # Generate plots if requested
    if plot_contacts and time_step is not None:
        plot_contact_heatmap(contact_matrix, skip_frames, time_step, x_labels, y_labels, atom_type,
                             output_dir=output_dir)
        plot_contact_distance(sub_dm_total, contact_matrix, time_step, skip_frames, output_dir=output_dir)

    if plot_distance_matrix:
        plot_distance_heatmap(sub_dm_total, x_labels, y_labels, atom_type, output_dir=output_dir)

    # Save data files if requested
    if save_data:
        save_matrix_data(sub_dm_total, contact_matrix, output_dir=output_dir)

    return sub_dm_total, contact_matrix
