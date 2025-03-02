# CRISP/data_analysis/cn_cal.py

from ase.io import read  
import numpy as np
from joblib import Parallel, delayed
import pickle
import argparse
from typing import Union, List, Dict, Tuple, Optional, Any
import os
import matplotlib.pyplot as plt
import itertools

def indices(atoms, ind: Union[str, List[Union[int, str]]]) -> np.ndarray:
    """
    Get atom indices based on different input types.
    
    Parameters:
        atoms: ASE Atoms object
            The atomic structure
        ind: Union[str, List[Union[int, str]]]
            Input specifier for atom indices. Can be:
            - "all" or None: all atoms
            - string ending with ".npy": path to numpy file with indices
            - integer or list of integers: directly as indices
            - string or list of strings: chemical symbols to select
    
    Returns:
        np.ndarray: Array of atom indices
        
    Raises:
        ValueError: If the input type is invalid
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
        for symbol in ind:
            idx.append(np.where(np.array(atoms.get_chemical_symbols()) == symbol)[0])
        return np.concatenate(idx)
    raise ValueError("Invalid index type")

def process_frame(frame, index_type) -> Tuple[np.ndarray, np.ndarray]:
    """
    Process a single frame to get distance matrix, with optimized handling for specific atom types.
    
    Parameters:
        frame: ASE Atoms object
            The atomic structure frame to process
        index_type: Union[str, List[Union[int, str]]]
            Specifier for which atom indices to include
    
    Returns:
        Tuple[np.ndarray, np.ndarray]: 
            - Complete distance matrix for the frame
            - Subset distance matrix filtered by the specified indices
    """
    dm = frame.get_all_distances(mic=True)
    idx = indices(frame, index_type)
    # Get sub-distance matrix (this will be used in other functions if needed)
    sub_dm = dm[np.ix_(idx, idx)]
    return dm, sub_dm

def calculate_coordination(frame, dm: np.ndarray, index_type, atom: str, cutoffs: Dict) -> Dict[int, int]:
    """
    Calculate coordination numbers with optimized handling for a single cutoff pair.
    
    Parameters:
        frame: ASE Atoms object
            The atomic structure frame
        dm: np.ndarray
            Distance matrix for the frame
        index_type: Union[str, List[Union[int, str]]]
            Specifier for which atom indices to include
        atom: str
            Chemical symbol of the atom type to calculate coordination for
        cutoffs: Dict[Tuple[str, str], float]
            Dictionary mapping a pair of atom symbols to cutoff distance
    
    Returns:
        Dict[int, int]: Dictionary mapping atom indices to their coordination numbers
        
    Raises:
        ValueError: If more than one cutoff pair is provided
    """
    # Verify that we only have a single cutoff pair
    if len(cutoffs) != 1:
        raise ValueError("This implementation only supports a single cutoff pair")
    
    # Extract the single cutoff value and symbols
    (symbol1, symbol2), cutoff_value = next(iter(cutoffs.items()))
    
    symbols = np.array(frame.get_chemical_symbols())
    
    # Get indices provided in index_type to restrict our search
    idx = indices(frame, index_type)
    
    # Only consider atoms in the provided indices
    idx_set = set(idx)  # For fast membership testing
    
    # Get target indices for specified atom type
    target_indices = [i for i in idx if symbols[i] == atom]
    
    coordination_dict = {}
    
    # If this is a homo-atomic pair (like O-O)
    homo_atomic = symbol1 == symbol2
    
    for i in target_indices:
        # Skip if atom type doesn't match either symbol in the cutoff pair
        if symbols[i] != symbol1 and symbols[i] != symbol2:
            coordination_dict[i] = 0
            continue
        
        # Find which atoms to count as neighbors (ONLY from the provided indices)
        if homo_atomic:
            # For homo-atomic pairs, count atoms of the same type but exclude self
            neighbor_mask = np.array([j in idx_set and symbols[j] == symbol1 and j != i for j in range(len(symbols))])
        else:
            # For hetero-atomic pairs, count atoms of the other type
            if symbols[i] == symbol1:
                neighbor_mask = np.array([j in idx_set and symbols[j] == symbol2 for j in range(len(symbols))])
            else:
                neighbor_mask = np.array([j in idx_set and symbols[j] == symbol1 for j in range(len(symbols))])
        
        # Get distances to eligible neighbor atoms
        distances = dm[i][neighbor_mask]
        
        # Count neighbors within cutoff distance
        coordination = np.sum(distances < cutoff_value)
        coordination_dict[i] = coordination
    
    return coordination_dict

def calculate_atom_coordination(
    file_path: str,
    output_dir: str,
    frame_skip: int,
    index_type: Union[str, List[Union[int, str]]],
    cutoffs: Dict[Tuple[str, str], float],
    atom: str,
    dm_path: Optional[str] = None,
    plot_data: bool = False) -> None:
    """
    Calculate coordination numbers for specified atom types across a trajectory.
    
    Parameters:
        file_path: str
            Path to the ASE trajectory file
        output_dir: str
            Directory to save output files
        frame_skip: int
            Number of frames to skip between calculations
        index_type: Union[str, List[Union[int, str]]]
            Specifier for which atom indices to include
        cutoffs: Dict[Tuple[str, str], float]
            Dictionary mapping a pair of atom symbols to cutoff distance
        atom: str
            Chemical symbol of the atom type to calculate coordination for
        dm_path: Optional[str]
            Path to precomputed distance matrices (if available)
        plot_data: bool, optional
            Whether to automatically generate plots (default: False)
    
    Returns:
        None: Results are saved to disk rather than returned
    """
    
    # Verify we only have a single cutoff pair
    if len(cutoffs) != 1:
        raise ValueError("This implementation only supports a single cutoff pair")
    
    # Use ase.io.read instead of Trajectory
    frames = read(file_path, index=f'::{frame_skip}')
    
    # Ensure frames is a list even for single frames
    if not isinstance(frames, list):
        frames = [frames]
        
    print(f"Processing trajectory with {len(frames)} frames (every {frame_skip}th frame)")
    
    # Create output directory if needed
    os.makedirs(output_dir, exist_ok=True)
    output_file = os.path.join(output_dir, "coordination_results.pkl")
    
    # Load or calculate distance matrices
    if dm_path:
        dms = np.load(dm_path, allow_pickle=True)
    else:
        results = Parallel(n_jobs=-1)(delayed(process_frame)(frame, index_type) for frame in frames)
        dms, _ = zip(*results)

    # Save distance matrices if we calculated them
    if not dm_path:
        dm_save_path = os.path.join(output_dir, "distance_matrices.npy")
        np.save(dm_save_path, dms)
        print(f"Saved distance matrices to {dm_save_path}")
    
    # Parallel coordination calculation
    coordination_types_list = Parallel(n_jobs=-1)(
        delayed(calculate_coordination)(frame, dm, index_type, atom, cutoffs) 
        for frame, dm in zip(frames, dms)
    )
    
    # Print the first frame's coordination data
    if coordination_types_list:
        print(f"Coordination data for first frame: {coordination_types_list[0]}")
    
    # Save results
    with open(output_file, 'wb') as f:
        pickle.dump(coordination_types_list, f)

    print(f"Analysis complete. Results saved to {output_file}!")

    # Generate and save plots if requested
    if plot_data:
        generate_coordination_plots(coordination_types_list, atom, output_dir)
        print(f"Plots generated and saved in {output_dir}")
    else:
        print("Plot generation skipped. To visualize results, use:")
        print(f"  - Load data: coordination_data = pickle.load(open('{output_file}', 'rb'))")
        print(f"  - Generate plots: generate_coordination_plots(coordination_data, '{atom}', '{output_dir}')")
    
    return

def parse_cutoffs(cutoff_str: str) -> Dict[Tuple[str, str], float]:
    """
    Parse a single cutoff pair from string format 'symbol1-symbol2:value'.
    
    Parameters:
        cutoff_str: str
            String specifying cutoff in format 'symbol1-symbol2:value'
    
    Returns:
        Dict[Tuple[str, str], float]: Dictionary mapping atom symbol pair to cutoff distance
        
    Raises:
        ValueError: If more than one cutoff pair is provided (contains commas)
    """
    cutoffs = {}
    if ',' in cutoff_str:
        raise ValueError("This implementation only supports a single cutoff pair")
    
    symbols_part, cutoff = cutoff_str.split(':')
    s1, s2 = symbols_part.split('-')
    cutoffs[(s1.strip(), s2.strip())] = float(cutoff)
    return cutoffs
                                     
                                     
def main() -> None:
    """
    Main function that handles command line arguments and runs coordination number calculation.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("trajectory", type=str)
    parser.add_argument("output_dir", type=str)
    parser.add_argument("--frame_skips", type=int, default=10)
    parser.add_argument("--index_type", nargs='+', default="all")
    parser.add_argument("--cutoffs", type=str, required=True)
    parser.add_argument("--atom", type=str, required=True)
    parser.add_argument("--dm_path", type=str)
    parser.add_argument("--generate_plots", action="store_true",  # New flag name that's more intuitive
                        help="Generate plots after analysis")

    args = parser.parse_args()

    # Create output directory if needed
    os.makedirs(args.output_dir, exist_ok=True)

    # Parse cutoffs
    cutoffs = parse_cutoffs(args.cutoffs)

    # Run coordination calculation
    calculate_atom_coordination(
        args.trajectory,
        args.output_dir,
        args.frame_skips,
        args.index_type,
        cutoffs,
        args.atom,
        args.dm_path,
        plot_data=args.generate_plots  # Use the flag directly now
    )


if __name__ == "__main__":
    main()

def calculate_avg_percentages(coordination_data):
    """
    Calculate the percentage of atoms with each coordination number across frames.
    
    Parameters:
        coordination_data: List[Dict[int, int]]
            List of dictionaries mapping atom indices to coordination numbers
    
    Returns:
        Dict[int, List[float]]: 
            Dictionary mapping coordination types to lists of percentage values across frames
    """
    # Get unique coordination types from the data
    coord_types = set(itertools.chain.from_iterable(frame_data.values() for frame_data in coordination_data))
    coord_types = sorted(coord_types)  # Sort the coordination types for consistent plotting
    avg_percentages = {coord_type: [] for coord_type in coord_types}

    for frame_data in coordination_data:
        total_atoms = len(frame_data)
        
        for coord_type in coord_types:
            count = sum(1 for value in frame_data.values() if value == coord_type)
            avg_percentage = count / total_atoms * 100
            avg_percentages[coord_type].append(avg_percentage)
    
    return avg_percentages

def plot_avg_percentages(avg_percentages, atom_type, output_file=None):
    """
    Plot the percentages of atoms with each coordination number across frames.
    
    Parameters:
        avg_percentages: Dict[int, List[float]]
            Dictionary mapping coordination types to lists of percentage values
        atom_type: str
            Chemical symbol of the atom type being analyzed
        output_file: str, optional
            Path to save the plot (if None, plot is displayed)
            
    Returns:
        None
    """
    frames = list(range(len(next(iter(avg_percentages.values())))))
    coord_types = list(avg_percentages.keys())
    
    # Generate colors and markers based on the number of coordination types
    colors = plt.cm.get_cmap('tab10', len(coord_types)).colors
    markers = itertools.cycle(['o', 's', 'D', '^', 'v', 'p', '*', '+', 'x'])

    plt.figure(figsize=(10, 6))

    for i, coord_type in enumerate(coord_types):
        plt.plot(frames, avg_percentages[coord_type], label=f'CN={coord_type}', 
                 color=colors[i], marker=next(markers), markevery=max(1, len(frames)//20))
                 
    # Calculate and plot mean lines
    for i, coord_type in enumerate(coord_types):
        mean_value = sum(avg_percentages[coord_type]) / len(avg_percentages[coord_type])
        plt.axhline(y=mean_value, color=colors[i], linestyle='--', alpha=0.7, 
                   label=f'Mean CN={coord_type}: {mean_value:.1f}%')

    plt.xlabel('Frame Index', fontsize=12)
    plt.ylabel('Percentage of Atoms (%)', fontsize=12)
    plt.title(f'Coordination Analysis: {atom_type} Atoms', fontsize=14)
    
    # Place legend outside the plot
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    
    if output_file:
        plt.savefig(output_file, bbox_inches='tight')  # Include legend outside
        plt.show()
        plt.close()

def calculate_overall_avg_percentages(avg_percentages):
    """
    Calculate the overall average percentages across all frames.
    
    Parameters:
        avg_percentages: Dict[int, List[float]]
            Dictionary mapping coordination types to lists of percentage values
    
    Returns:
        Dict[int, float]: Dictionary mapping coordination types to overall average percentages
    """
    return {coord_type: sum(percentages) / len(percentages) for coord_type, percentages in avg_percentages.items()}

def plot_coordination_distribution(coordination_data, atom_type, output_file=None):
    """
    Create a pie chart showing the average distribution of coordination numbers.
    
    Parameters:
        coordination_data: List[Dict[int, int]]
            List of dictionaries mapping atom indices to coordination numbers
        atom_type: str
            Chemical symbol of the atom type being analyzed
        output_file: str, optional
            Path to save the plot
            
    Returns:
        None
    """
    avg_percentages = calculate_avg_percentages(coordination_data)
    overall_avg_percentages = calculate_overall_avg_percentages(avg_percentages)
    
    # Create a pie chart
    plt.figure(figsize=(10, 7))
    
    labels = [f'CN={coord_type}: {pct:.1f}%' for coord_type, pct in overall_avg_percentages.items()]
    plt.pie(overall_avg_percentages.values(), 
            labels=labels, 
            autopct='%1.1f%%',
            colors=plt.cm.tab10.colors[:len(overall_avg_percentages)],
            startangle=90)
    
    plt.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle
    plt.title(f'Average Distribution of {atom_type} Atom Coordination', fontsize=14)
    
    if output_file:
        plt.savefig(output_file, bbox_inches='tight')
        plt.show()
        plt.close()

def generate_coordination_plots(coordination_data, atom_type, output_dir):
    """
    Generate and save plots for coordination analysis.
    
    Parameters:
        coordination_data: List[Dict[int, int]]
            List of dictionaries mapping atom indices to coordination numbers
        atom_type: str
            Chemical symbol of the atom type being analyzed
        output_dir: str
            Directory to save the plots
            
    Returns:
        None
    """
    # Calculate percentages for each coordination number across frames
    avg_percentages = calculate_avg_percentages(coordination_data)
    
    # Create time series plot
    time_series_file = os.path.join(output_dir, f"{atom_type}_coordination_time_series.png")
    plot_avg_percentages(avg_percentages, atom_type, time_series_file)
    
    # Create distribution pie chart
    pie_chart_file = os.path.join(output_dir, f"{atom_type}_coordination_distribution.png")
    plot_coordination_distribution(coordination_data, atom_type, pie_chart_file)
    
    # Calculate and save overall statistics
    overall_avg_percentages = calculate_overall_avg_percentages(avg_percentages)
    
    # Save statistics to text file
    stats_file = os.path.join(output_dir, f"{atom_type}_coordination_statistics.txt")
    with open(stats_file, 'w') as f:
        f.write(f"Coordination Analysis for {atom_type} Atoms\n")
        f.write("======================================\n\n")
        f.write("Overall Average Percentages:\n")
        for coord_type, avg_percentage in sorted(overall_avg_percentages.items()):
            f.write(f"  CN={coord_type}: {avg_percentage:.2f}%\n")
            
    print(f"Generated coordination plots for {atom_type} atoms:")
    print(f"  - Time series plot: {time_series_file}")
    print(f"  - Distribution pie chart: {pie_chart_file}")
    print(f"  - Statistics file: {stats_file}")