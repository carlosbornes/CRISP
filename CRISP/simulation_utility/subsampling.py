# CRISP/simulation_utility/subsampling.py

import numpy as np
from ase.io import read, write
import fpsample
import glob
from dscribe.descriptors import SOAP
import matplotlib.pyplot as plt
import argparse
import os

def indices(atoms, ind):
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
    numpy.ndarray
        Array of atom indices
    """
    if ind == "all" or ind is None:
        idx = np.arange(len(atoms))
        return idx  # Select all atoms by default

    if isinstance(ind, str) and ind.endswith(".npy"):
        idx = np.load(ind, allow_pickle=True)
        return idx  # For numpy file with indices

    if not isinstance(ind, list): 
        ind = [ind]

    if any(isinstance(item, int) for item in ind):
        idx = np.array(ind)
        return idx  # For specific atoms like [1,2,3] or 1

    if any(isinstance(item, str) for item in ind):
        idx = []
        for i in ind:
            idx.append(np.where(np.array(atoms.get_chemical_symbols()) == i)[0])
        idx = np.concatenate(idx)
        return idx  # For chemical symbols - either "O", ["O", "H"] or ("O", "H")

def create_repres(traj, ind="all"):
    """Create SOAP representation vectors for trajectory frames.
    
    Parameters
    ----------
    traj : list
        List of ASE Atoms objects (trajectory frames)
    ind : str, list, or None, optional
        Specification for which atoms to select for SOAP calculation (default: "all")
        
    Returns
    -------
    numpy.ndarray
        Array of SOAP feature vectors, one per frame
    """
    all_spec = traj[0].get_chemical_symbols()
    idx = indices(traj[0], ind=ind)
    
    repres = []
    for index, structure in enumerate(traj):
        periodic_cell = structure.cell.volume > 0

        soap = SOAP(
            species=all_spec,
            periodic=periodic_cell,
            r_cut=6,
            n_max=8,
            l_max=6,
            sigma=0.5,
            sparse=False
        )

        soap_ind = soap.create(traj[index], centers=idx)
        repres.append(np.mean(soap_ind, axis=0))

    return np.array(repres)

def run_subsample(filename, output_dir, n_samples=50, index_type="all", file_format=None, skip=1):
    """Perform FPS (Farthest Point Sampling) on a trajectory.
    
    Parameters
    ----------
    filename : str
        Path to trajectory file (can include wildcards)
    output_dir : str
        Directory to save output files
    n_samples : int, optional
        Number of structures to subsample (default: 50)
    index_type : str, list, or None, optional
        Atom indices to focus the subsampling on (default: "all")
    file_format : str, optional
        File format of the trajectories (default: None, auto-detect)
    skip : int, optional
        Read every n-th frame from trajectory (default: 1)
        
    Returns
    -------
    list
        List of selected ASE Atoms objects (subsampled frames)
    """
    os.makedirs(output_dir, exist_ok=True)  # Ensure output directory exists
    filename = glob.glob(filename)

    trajec = []
    for file in filename:
        if file_format is not None:
            trajec = read(file, index=f'::{skip}', format=file_format)
        else:
            trajec = read(file, index=f'::{skip}')

    if not isinstance(trajec, list): 
        trajec = [trajec]

    repres = create_repres(trajec, ind=index_type)
    perm = fpsample.fps_sampling(repres, n_samples, start_idx=0)

    # Create output trajectory file path
    traj_output_file = os.path.join(output_dir, 'subsampled_trajectory.xyz')
    
    # Delete the file if it already exists to avoid appending to existing file
    if os.path.exists(traj_output_file):
        os.remove(traj_output_file)
    
    fps_frames = []
    for str_idx, frame in enumerate(perm):
        new_frame = trajec[frame]
        fps_frames.append(new_frame)
        
        # Write to trajectory file, append=True after the first frame
        write(traj_output_file, new_frame, append=(str_idx > 0))
    
    print(f"Successfully wrote {n_samples} frames to {traj_output_file}")
    return fps_frames

def main():
    """Parse command-line arguments and run subsampling."""
    parser = argparse.ArgumentParser()
    parser.add_argument('filename', type=str)
    parser.add_argument('output_dir', type=str)
    parser.add_argument('--n_samples', type=int, default=50)
    parser.add_argument('--index_type', type=str, default='all')
    parser.add_argument('--file_format', type=str, default=None)
    parser.add_argument('--skip', type=int, default=1)
    
    args = parser.parse_args()

    run_subsample(args.filename, args.output_dir, args.n_samples, args.index_type, args.file_format, args.skip)

if __name__ == "__main__":
    main()
