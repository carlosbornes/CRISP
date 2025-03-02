# CRISP/data_analysis/h_bond.py

from ase.io.trajectory import Trajectory
import numpy as np
import csv
from joblib import Parallel, delayed
import argparse
import os
from typing import Union, List, Optional, Tuple, Any

def indices(atoms, ind=None, default_element=None):
    """
    Extract atom indices based on various input types.
    
    Parameters
    ----------
    atoms : ase.Atoms
        Atoms object containing atomic coordinates and elements
    ind : str, list, or None, optional
        Specification for which atoms to select:
        - None: use default_element to select atoms
        - string ending with ".npy": load indices from NumPy file
        - list of integers: direct atom indices
        - list of strings: chemical symbols to select
    default_element : str, optional
        Chemical symbol to use when ind is None
        
    Returns
    -------
    np.ndarray
        Array of atom indices
        
    Raises
    ------
    ValueError
        If the index type is not recognized
    """
    if ind is None and default_element is not None:
        return np.where(np.array(atoms.get_chemical_symbols()) == default_element)[0]
        
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

def get_distance_matrices(atoms, o_ind, h_ind):
    """
    Calculate distance matrices for oxygen-oxygen and hydrogen-oxygen pairs.
    
    Parameters
    ----------
    atoms : ase.Atoms
        Atoms object containing atomic coordinates
    o_ind : np.ndarray
        Indices of oxygen atoms
    h_ind : np.ndarray
        Indices of hydrogen atoms
        
    Returns
    -------
    Tuple[np.ndarray, np.ndarray, np.ndarray]
        Tuple containing:
        - full_dm: Full distance matrix for all atoms
        - o_o_dm: Distance matrix between oxygen atoms
        - h_o_dm: Distance matrix between hydrogen and oxygen atoms
    """
    full_dm = atoms.get_all_distances(mic=True)
    
    o_o_dm = full_dm[np.ix_(o_ind, o_ind)]
    
    h_o_dm = full_dm[np.ix_(h_ind, o_ind)]
    
    return full_dm, o_o_dm, h_o_dm

def count_hydrogen_bonds(atoms, o_ind, h_ind, oxygen_oxygen_distance_threshold, 
                        oxygen_hydrogen_distance_threshold, angle_cutoff, distance_matrices=None):
    """
    Count hydrogen bonds using distance matrices for efficiency.
    
    Parameters
    ----------
    atoms : ase.Atoms
        Atoms object containing atomic coordinates
    o_ind : np.ndarray
        Indices of oxygen atoms
    h_ind : np.ndarray
        Indices of hydrogen atoms
    oxygen_oxygen_distance_threshold : float
        Maximum distance between donor and acceptor oxygen atoms (Å)
    oxygen_hydrogen_distance_threshold : float
        Maximum distance between oxygen and hydrogen atoms (Å)
    angle_cutoff : float
        Maximum angle in degrees for hydrogen bond formation
    distance_matrices : tuple, optional
        If provided, (full_dm, o_o_dm, h_o_dm) matrices to use
        
    Returns
    -------
    Tuple[List[List], int]
        Tuple containing:
        - List of hydrogen bond data [h, o_donor, o_acceptor, r_oo, r_ho, angle]
        - Number of hydrogen bonds found
    """
    data = []
    num_hydrogen_bonds = 0
    used_hydrogens = set()  # Track hydrogens already involved in hydrogen bonds
    
    # Get distance matrices if not provided
    if distance_matrices is None:
        full_dm, o_o_dm, h_o_dm = get_distance_matrices(atoms, o_ind, h_ind)
    else:
        full_dm, o_o_dm, h_o_dm = distance_matrices
    
    # Build a map of donor oxygens to their attached hydrogens using the distance matrix
    donor_hydrogens = {}
    for i, o in enumerate(o_ind):
        donor_hydrogens[o] = []
        for j, h in enumerate(h_ind):
            if h_o_dm[j, i] < oxygen_hydrogen_distance_threshold:
                donor_hydrogens[o].append(h)

    # Check each potential donor-acceptor pair
    for i, o_donor in enumerate(o_ind):
        hydrogens = donor_hydrogens.get(o_donor, [])
        for h in hydrogens:
            if h in used_hydrogens:
                continue  # Skip hydrogens already in hydrogen bonds
            
            for j, o_acceptor in enumerate(o_ind):
                if o_acceptor == o_donor:
                    continue  # Skip the same oxygen
                
                r_oo = o_o_dm[i, j]  # Use pre-computed distance
                
                if r_oo < oxygen_oxygen_distance_threshold:
                    # Calculate angle using atoms object - this requires atom positions
                    angle = atoms.get_angle(o_acceptor, o_donor, h, mic=True)
                    
                    if angle < angle_cutoff:
                        # Hydrogen bond detected
                        r_ho = full_dm[h, o_acceptor]  # Use pre-computed distance
                        data.append([h, o_donor, o_acceptor, r_oo, r_ho, angle])
                        num_hydrogen_bonds += 1
                        used_hydrogens.add(h)  # Mark hydrogen as used
                        break  # A hydrogen can only form one H-bond

    return data, num_hydrogen_bonds

def analyze_hydrogen_bonds(
    traj_path: str,
    output_dir: str,
    oxygen_indices: Optional[Union[str, List[int], np.ndarray]] = None,
    hydrogen_indices: Optional[Union[str, List[int], np.ndarray]] = None,
    frame_skip: int = 1,
    oxygen_oxygen_distance_threshold: float = 3.5,
    oxygen_hydrogen_distance_threshold: float = 1.2,
    angle_cutoff: float = 30.0,
    dm_path: Optional[str] = None
):
    """
    Analyze hydrogen bonds in a molecular dynamics trajectory.
    
    Parameters
    ----------
    traj_path : str
        Path to the ASE trajectory file
    output_dir : str
        Directory to save output files
    oxygen_indices : Optional[Union[str, List[int], np.ndarray]], optional
        Oxygen atom indices (default: None, uses all O atoms)
    hydrogen_indices : Optional[Union[str, List[int], np.ndarray]], optional
        Hydrogen atom indices (default: None, uses all H atoms)
    frame_skip : int, optional
        Number of frames to skip between analyses (default: 1)
    oxygen_oxygen_distance_threshold : float, optional
        Maximum O-O distance for hydrogen bonds (default: 3.5 Å)
    oxygen_hydrogen_distance_threshold : float, optional
        Maximum O-H distance for hydrogen bonds (default: 1.2 Å)
    angle_cutoff : float, optional
        Maximum angle cutoff in degrees (default: 30.0°)
    dm_path : Optional[str], optional
        Path to pre-computed distance matrices (default: None)
        
    Returns
    -------
    List[int]
        Total number of hydrogen bonds per analyzed frame
    """
    os.makedirs(output_dir, exist_ok=True)
    
    # Set output filenames
    output_filename = os.path.join(output_dir, f'hydrogen_bonds_{frame_skip}skips.csv')
    total_bonds_filename = os.path.join(output_dir, f'total_hydrogen_bonds_per_frame_{frame_skip}skips.csv')
    
    # Load trajectory
    traj = Trajectory(traj_path)
    
    # Get frames to analyze
    frames = traj[::frame_skip]
    print(f"Analyzing {len(frames)} frames from trajectory")
    
    # Get first frame to determine indices
    first_frame = traj[0]
    
    # Get oxygen and hydrogen indices using the flexible indices function
    oxygen_idx = indices(first_frame, oxygen_indices, default_element='O')
    hydrogen_idx = indices(first_frame, hydrogen_indices, default_element='H')
    
    print(f"Using {len(oxygen_idx)} oxygen atoms for hydrogen bond analysis")
    print(f"Using {len(hydrogen_idx)} hydrogen atoms for hydrogen bond analysis")
    
    # Load or calculate distance matrices
    distance_matrices = None
    if dm_path:
        try:
            distance_matrices = np.load(dm_path, allow_pickle=True)
            print(f"Loaded distance matrices from {dm_path}")
        except Exception as e:
            print(f"Error loading distance matrices: {e}")
            print("Will calculate distance matrices instead")
    
    # Parallelize the hydrogen bond analysis using joblib
    all_data = Parallel(n_jobs=-1)(delayed(count_hydrogen_bonds)(
        atoms, oxygen_idx, hydrogen_idx,
        oxygen_oxygen_distance_threshold, oxygen_hydrogen_distance_threshold, angle_cutoff,
        distance_matrices[i] if distance_matrices is not None else None
    ) for i, atoms in enumerate(frames))
    
    # Extract data and total hydrogen bonds per frame
    total_hydrogen_bonds_per_frame = [num_bonds for _, num_bonds in all_data]
    
    # Create data with frame indices
    data_with_frames = []
    for frame_idx, (sublist, _) in enumerate(all_data):
        actual_frame_idx = frame_idx * frame_skip
        for entry in sublist:
            data_with_frames.append([actual_frame_idx] + entry)
    
    # Save the hydrogen bond data to a CSV file with frame information
    with open(output_filename, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['Frame', 'Hydrogen', 'Donor', 'Acceptor', 
                        'Oxygen-Oxygen Distance', 'Hydrogen-Oxygen Distance', 'Angle(OA-OD-HD)'])
        writer.writerows(data_with_frames)
    
    # Save the total hydrogen bonds per frame to a separate CSV file
    with open(total_bonds_filename, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['Frame', 'Total Hydrogen Bonds'])
        for frame_idx, num_bonds in enumerate(total_hydrogen_bonds_per_frame):
            actual_frame_idx = frame_idx * frame_skip
            writer.writerow([actual_frame_idx, num_bonds])
    
    print(f"Analysis complete. Results saved to {output_dir}")
    return total_hydrogen_bonds_per_frame

def main():
    """
    Parse command-line arguments and run hydrogen bond analysis.
    
    Command-line Parameters
    ----------------------
    trajectory : str
        Path to ASE trajectory file
    output_dir : str
        Directory to save output files
    --oxygen_indices : str, optional
        Oxygen atoms to analyze
    --hydrogen_indices : str, optional
        Hydrogen atoms to analyze
    --frame_skip : int, optional
        Number of frames to skip between analyses
    --oo_threshold : float, optional
        Oxygen-oxygen distance threshold in Å
    --oh_threshold : float, optional
        Oxygen-hydrogen distance threshold in Å
    --angle_cutoff : float, optional
        Maximum angle cutoff in degrees
    --dm_path : str, optional
        Path to pre-computed distance matrices
        
    Returns
    -------
    None
    """
    parser = argparse.ArgumentParser()
    
    # Required arguments
    parser.add_argument("trajectory", type=str)
    parser.add_argument("output_dir", type=str)
    
    # Atom selection options
    atom_group = parser.add_argument_group('Atom selection')
    atom_group.add_argument("--oxygen_indices", type=str)
    atom_group.add_argument("--hydrogen_indices", type=str)
    
    # Analysis parameters
    analysis_group = parser.add_argument_group('Analysis parameters')
    analysis_group.add_argument("--frame_skip", type=int, default=100)
    analysis_group.add_argument("--oo_threshold", type=float, default=3.5)
    analysis_group.add_argument("--oh_threshold", type=float, default=1.2)
    analysis_group.add_argument("--angle_cutoff", type=float, default=30.0)
    analysis_group.add_argument("--dm_path", type=str)
    
    args = parser.parse_args()
    
    # Validate trajectory file
    if not os.path.exists(args.trajectory):
        raise FileNotFoundError(f"Trajectory file not found: {args.trajectory}")
    
    # Handle potential input conversions
    oxygen_indices = args.oxygen_indices
    hydrogen_indices = args.hydrogen_indices
    
    # Try to convert string representations of integers or lists
    if oxygen_indices:
        # Check if it's a file first
        if oxygen_indices.endswith('.npy') and os.path.exists(oxygen_indices):
            pass  # Keep as is - it's a valid file path
        # Try to interpret as list of integers
        elif oxygen_indices.startswith('[') and oxygen_indices.endswith(']'):
            try:
                oxygen_indices = eval(oxygen_indices)
            except:
                # If not a valid list format, treat as a symbol
                if '[' not in oxygen_indices[1:-1] and ',' not in oxygen_indices[1:-1]:
                    oxygen_indices = oxygen_indices.strip('[]')
            
    if hydrogen_indices:
        # Check if it's a file first
        if hydrogen_indices.endswith('.npy') and os.path.exists(hydrogen_indices):
            pass  # Keep as is - it's a valid file path
        # Try to interpret as list of integers
        elif hydrogen_indices.startswith('[') and hydrogen_indices.endswith(']'):
            try:
                hydrogen_indices = eval(hydrogen_indices)
            except:
                # If not a valid list format, treat as a symbol
                if '[' not in hydrogen_indices[1:-1] and ',' not in hydrogen_indices[1:-1]:
                    hydrogen_indices = hydrogen_indices.strip('[]')
    
    # Run the analysis
    analyze_hydrogen_bonds(
        args.trajectory,
        args.output_dir,
        oxygen_indices=oxygen_indices,
        hydrogen_indices=hydrogen_indices,
        frame_skip=args.frame_skip,
        oxygen_oxygen_distance_threshold=args.oo_threshold,
        oxygen_hydrogen_distance_threshold=args.oh_threshold,
        angle_cutoff=args.angle_cutoff,
        dm_path=args.dm_path
    )

if __name__ == "__main__":
    main()
