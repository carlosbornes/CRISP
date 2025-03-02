# CRISP/simulation_utility/atomic_indices.py

import os
import numpy as np
import ase.io
import csv    
import argparse


def atom_indices(file_path, frame_index=0, cutoffs=None):
    """Extract atom indices by chemical symbol and find atom pairs within specified cutoffs.
    
    Parameters
    ----------
    file_path : str
        Path to the ASE trajectory file
    frame_index : int, optional
        Index of the frame to analyze (default: 0)
    cutoffs : dict, optional
        Dictionary with atom symbol pairs as keys and cutoff distances as values
        Example: {('Si', 'O'): 2.0, ('Al', 'O'): 2.1}
        
    Returns
    -------
    indices_by_symbol : dict
        Dictionary with chemical symbols as keys and lists of atomic indices as values
    dist_matrix : numpy.ndarray
        Distance matrix between all atoms, accounting for periodic boundary conditions
    cutoff_indices : dict
        Dictionary with atom symbol pairs as keys and lists of (idx1, idx2, distance) tuples
        for atoms that are within the specified cutoff distance
    """
    new = ase.io.read(file_path, index=frame_index, format="traj")
    dist_matrix = new.get_all_distances(mic=True)
    symbols = new.get_chemical_symbols()

    unique_symbols = list(set(symbols))

    indices_by_symbol = {symbol: [] for symbol in unique_symbols}

    for idx, atom in enumerate(new):
        indices_by_symbol[atom.symbol].append(idx)

    cutoff_indices = {}

    if cutoffs:
        for pair, cutoff in cutoffs.items():
            symbol1, symbol2 = pair
            pair_indices_distances = []
            if symbol1 in indices_by_symbol and symbol2 in indices_by_symbol:
                for idx1 in indices_by_symbol[symbol1]:
                    for idx2 in indices_by_symbol[symbol2]:
                        if dist_matrix[idx1, idx2] < cutoff:
                            pair_indices_distances.append(
                                (idx1, idx2, dist_matrix[idx1, idx2])
                            )
            cutoff_indices[pair] = pair_indices_distances

    return indices_by_symbol, dist_matrix, cutoff_indices


def run_atom_indices(file_path, output_folder, frame_index=0, cutoffs=None):
    """Run atom index extraction and save results to files.
    
    Parameters
    ----------
    file_path : str
        Path to the ASE trajectory file
    output_folder : str
        Directory where output files will be saved
    frame_index : int, optional
        Index of the frame to analyze (default: 0)
    cutoffs : dict, optional
        Dictionary with atom symbol pairs as keys and cutoff distances as values
        
    Returns
    -------
    None
        Results are saved to the specified output folder:
        - lengths.npy: Dictionary of number of atoms per element
        - {symbol}_indices.npy: Numpy array of atom indices for each element
        - cutoff/{symbol1}-{symbol2}_cutoff.csv: CSV files with atom pairs within cutoff
    """
    indices, dist_matrix, cutoff_indices = atom_indices(file_path, frame_index, cutoffs)

    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    lengths = {symbol: len(indices[symbol]) for symbol in indices}
    np.save(os.path.join(output_folder, "lengths.npy"), lengths)

    for symbol, data in indices.items():
        np.save(os.path.join(output_folder, f"{symbol}_indices.npy"), data)
        print(f"Length of {symbol} indices: {len(data)}")

    print("Outputs saved.")

    # Save cutoff indices into CSV files
    cutoff_folder = os.path.join(output_folder, "cutoff")
    if not os.path.exists(cutoff_folder):
        os.makedirs(cutoff_folder)

    for pair, pair_indices_distances in cutoff_indices.items():
        symbol1, symbol2 = pair
        filename = f"{symbol1}-{symbol2}_cutoff.csv"
        filepath = os.path.join(cutoff_folder, filename)
        with open(filepath, mode="w", newline="") as file:
            writer = csv.writer(file)
            writer.writerow([f"{symbol1} index", f"{symbol2} index", "distance"])
            writer.writerows(pair_indices_distances)
        print(f"Saved cutoff indices for {symbol1}-{symbol2} to {filepath}")
        

def main():
    """Parse command-line arguments and run atom indices extraction."""
    parser = argparse.ArgumentParser()
    parser.add_argument("file_path", type=str)
    parser.add_argument("output_folder", type=str)
    parser.add_argument(
        "--frame_index",
        type=int,
        default=0,
    )
    parser.add_argument(
        "--cutoffs",
        type=str,
        default=None,
    )

    args = parser.parse_args()
    
    # Validate frame index
    try:
        traj_length = len(ase.io.read(args.file_path, index=":", format="traj"))
        
        # Check if frame_index is within valid range
        if args.frame_index < 0 or args.frame_index >= traj_length:
            raise ValueError(
                f"Error: Frame index {args.frame_index} is out of range. "
                f"Valid range is 0 to {traj_length-1}."
            )
            
        frame_index = args.frame_index
        print(f"Analyzing frame with index {frame_index} (out of {traj_length} frames)")
        
    except ValueError as e:
        # Re-raise ValueError for frame index out of range
        raise e
        
    except Exception as e:
        # Handle other exceptions when reading the trajectory
        raise ValueError(f"Error reading trajectory: {e}")

    cutoffs = None
    if args.cutoffs:
        cutoffs = {}
        pairs = args.cutoffs.split(",")
        for pair in pairs:
            symbols, cutoff = pair.split(":")
            symbol1, symbol2 = symbols.split("-")
            cutoffs[(symbol1, symbol2)] = float(cutoff)

    run_atom_indices(args.file_path, args.output_folder, frame_index, cutoffs)


if __name__ == "__main__":
    main()
