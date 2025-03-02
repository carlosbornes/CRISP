# CRISP/simulation_utility/interatomic_distances.py

import os
import numpy as np
import pickle
import argparse
from ase.io import read
from ase import Atoms
from typing import Union, Tuple, List, Dict, Any
from joblib import Parallel, delayed


def indices(atoms: Atoms, ind: Union[str, List[Union[int, str]], None]) -> np.ndarray:
    """Extract atom indices from various input types.
    
    Parameters
    ----------
    atoms : ase.Atoms
        Atoms object containing atomic coordinates and elements
    ind : str, list, or None
        Specification for which atoms to select:
        - "all" or None: all atoms
        - string ending with ".npy": load indices from NumPy file
        - list of integers: direct atom indices
        - list of strings: chemical symbols to select
        
    Returns
    -------
    np.ndarray
        Array of atom indices
        
    Raises
    ------
    ValueError
        If the index type is not recognized
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


def distance_calculation(
    traj_path: str,
    frame_skip: int,
    index_type: Union[str, List[Union[int, str]]] = "all"
) -> Tuple[List[np.ndarray], List[np.ndarray]]:
    """Calculate distance matrices for multiple frames in a trajectory.
    
    Parameters
    ----------
    traj_path : str
        Path to the trajectory file
    frame_skip : int
        Read every nth frame (n=frame_skip)
    index_type : str, list, or None, optional
        Specification for which atoms to select for sub-matrix (default: "all")
        
    Returns
    -------
    Tuple[List[np.ndarray], List[np.ndarray]]
        Two lists containing:
        1. Full distance matrices for all frames
        2. Sub-matrices for specified atoms
        
    Raises
    ------
    ValueError
        If no frames were found in the trajectory
    """
    frames = read(traj_path, index=f"::{frame_skip}")
    if not frames:
        raise ValueError("No frames were found in the trajectory using the given frame_skip.")

    def process_frame(frame: Atoms) -> Tuple[np.ndarray, np.ndarray]:
        dm = frame.get_all_distances(mic=True)
        idx = indices(frame, index_type)
        sub_dm = dm[np.ix_(idx, idx)]
        return dm, sub_dm

    results = Parallel(n_jobs=-1)(delayed(process_frame)(frame) for frame in frames)
    full_dms, sub_dms = zip(*results)
    return list(full_dms), list(sub_dms)


def save_distance_matrices(
    full_dms: List[np.ndarray], 
    sub_dms: List[np.ndarray],
    index_type: Union[str, List[Union[int, str]]] = "all",
    output_dir: str = "distance_calculations"
) -> None:
    """Save distance matrices to pickle file.
    
    Parameters
    ----------
    full_dms : List[np.ndarray]
        List of full distance matrices
    sub_dms : List[np.ndarray]
        List of sub-matrices for specified atoms
    index_type : str, list, or None, optional
        Type of index selection used (default: "all")
    output_dir : str, optional
        Directory to save output file (default: "distance_calculations")
        
    Returns
    -------
    None
        Saves results to disk
    """
    data = {"full_dms": full_dms}
    if index_type not in ["all", None]:
        data["sub_dms"] = sub_dms
        
    os.makedirs(output_dir, exist_ok=True)
    output_path = os.path.join(output_dir, "distance_matrices.pkl")
    with open(output_path, "wb") as f:
        pickle.dump(data, f)
    print(f"Distance matrices saved in '{output_path}'")


def main():
    """Parse command-line arguments and run distance calculations."""
    parser = argparse.ArgumentParser()
    parser.add_argument("traj_path", type=str)
    
    parser.add_argument(
        "--frame_skip",
        type=int,
        default=10,
    )
    
    parser.add_argument(
        "--index_type",
        type=str,
        default=None
    )
    
    parser.add_argument(
        "--output_dir", 
        type=str,
        default="distance_calculations"
    )

    args = parser.parse_args()
    
    full_dms, sub_dms = distance_calculation(args.traj_path, args.frame_skip, args.index_type) 
    save_distance_matrices(full_dms, sub_dms, args.index_type, args.output_dir)


if __name__ == "__main__":
    main()
