"""
CRISP/simulation_utility/subsampling.py

This module provides functionality for structure subsampling from molecular dynamics
trajectories using Farthest Point Sampling (FPS) with SOAP descriptors.
"""

import numpy as np
from ase.io import read, write
import fpsample
import glob
import os
from dscribe.descriptors import SOAP
import matplotlib.pyplot as plt
from joblib import Parallel, delayed


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
        if isinstance(ind, str):
            ind = [ind]
        for i in ind:
            idx.append(np.where(np.array(atoms.get_chemical_symbols()) == i)[0])
        idx = np.concatenate(idx)
        return idx  # For chemical symbols - either "O", ["O", "H"] or ("O", "H")


def compute_soap(structure, all_spec, rcut, idx):
    """Compute SOAP descriptors for a given structure.
    
    Parameters
    ----------
    structure : ase.Atoms
        Atomic structure for which to compute SOAP descriptors
    all_spec : list
        List of chemical elements to include in the descriptor
    rcut : float
        Cutoff radius for the SOAP descriptor in Angstroms
    idx : numpy.ndarray
        Indices of atoms to use as centers for SOAP calculation
        
    Returns
    -------
    numpy.ndarray
        Average SOAP descriptor vector for the structure
    """
    periodic_cell = structure.cell.volume > 0
    soap = SOAP(
        species=all_spec,
        periodic=periodic_cell,
        r_cut=rcut,
        n_max=8,
        l_max=6,
        sigma=0.5,
        sparse=False
    )
    soap_ind = soap.create(structure, centers=idx)
    return np.mean(soap_ind, axis=0)


def create_repres(traj, rcut=6, ind="all", n_jobs=-1):
    """Create SOAP representation vectors for a trajectory.
    
    Parameters
    ----------
    traj : list
        List of ase.Atoms objects representing a trajectory
    rcut : float, optional
        Cutoff radius for the SOAP descriptor in Angstroms (default: 6)
    ind : str, list, or None, optional
        Specification for which atoms to use as SOAP centers (default: "all")
    n_jobs : int, optional
        Number of parallel jobs to run; -1 uses all available cores (default: -1)
        
    Returns
    -------
    numpy.ndarray
        Array of SOAP descriptors for each frame in the trajectory
    """
    all_spec = traj[0].get_chemical_symbols()
    idx = indices(traj[0], ind=ind)

    repres = Parallel(n_jobs=n_jobs)(
        delayed(compute_soap)(structure, all_spec, rcut, idx) for structure in traj
    )

    return np.array(repres)


def subsample(filename, n_samples=50, index_type="all", rcut=6.0, file_format=None, 
             plot_subsample=False, skip=1, output_dir="subsampled_structures"):
    """Subsample a trajectory using Farthest Point Sampling with SOAP descriptors.
    
    Parameters
    ----------
    filename : str
        Path pattern to trajectory file(s); supports globbing
    n_samples : int, optional
        Number of frames to select (default: 50)
    index_type : str, list, or None, optional
        Specification for which atoms to use for SOAP calculation (default: "all")
    rcut : float, optional
        Cutoff radius for SOAP in Angstroms (default: 6.0)
    file_format : str, optional
        File format for ASE I/O (default: None, auto-detect)
    plot_subsample : bool, optional
        Whether to generate a plot of FPS distances (default: False)
    skip : int, optional
        Read every nth frame from the trajectory (default: 1)
    output_dir : str, optional
        Directory to save the subsampled structures (default: "subsampled_structures")
        
    Returns
    -------
    list
        List of selected ase.Atoms frames
        
    Notes
    -----
    The selected frames and plots are saved in the specified output directory
    """
    filename = glob.glob(filename)

    trajec = []
    for file in filename:
        if file_format is not None:
            trajec += read(file, index=f'::{skip}', format=file_format)
        else:
            trajec += read(file, index=f'::{skip}')

    if not isinstance(trajec, list): 
        trajec = [trajec]
    
    repres = create_repres(trajec, ind=index_type, rcut=rcut)
    perm = fpsample.fps_sampling(repres, n_samples, start_idx=0)

    fps_frames = []

    for str_idx, frame in enumerate(perm):
        new_frame = trajec[frame]
        fps_frames.append(new_frame)

    # Create output directory if it doesn't exist
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    if plot_subsample:
        distance = []
        for i in range(1, len(perm)):
            distance.append(np.min(np.linalg.norm(repres[perm[:i]] - repres[perm[i]], axis=1)))

        plt.figure(figsize=(8, 6))
        plt.plot(distance, c="blue", linewidth=2)
        plt.ylim([0, 1.1 * max(distance)])
        plt.xlabel("Number of subsampled structures")
        plt.ylabel("Euclidean distance")
        plt.title("FPS Subsampling")
        # Save plot in the output directory
        plt.savefig(os.path.join(output_dir, "subsampled_convergence.png"), dpi=300)
        plt.show()
        plt.close()
        print(f"Saved convergence plot to {os.path.join(output_dir, 'subsampled_convergence.png')}")

    # Extract the base filename without path for output file
    base_filename = filename[0].split('/')[-1] if '/' in filename[0] else filename[0].split('\\')[-1]
    output_file = os.path.join(output_dir, f"subsample_{base_filename}")
    write(output_file, fps_frames, format=file_format)
    print(f"Saved {len(fps_frames)} subsampled structures to {output_file}")

    return fps_frames
