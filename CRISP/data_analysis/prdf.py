import os
import numpy as np
import pickle
import math
import argparse
from ase.io import read
from ase import Atoms
from typing import Optional, Union, Tuple, List
from joblib import Parallel, delayed


def check_cell_and_r_max(atoms: Atoms, rmax: float):
    """
    Check that the cell is large enough to contain a sphere of radius rmax.
    
    Parameters
    ----------
    atoms : Atoms
        ASE Atoms object with cell information
    rmax : float
        Maximum radius to consider
        
    Raises
    ------
    VolumeNotDefined
        If cell is not defined or too small for requested rmax
    """
    if not atoms.cell.any():
        raise VolumeNotDefined("The system's cell is not defined.")
    cell = atoms.cell
    try:
        lengths = cell.lengths()
        if np.min(lengths) < 2 * rmax:
            raise VolumeNotDefined(f"Cell length {np.min(lengths)} is smaller than 2*rmax ({2*rmax}).")
    except AttributeError:
        volume = cell.volume
        required_volume = (4/3) * math.pi * rmax**3
        if volume < required_volume:
            raise VolumeNotDefined(f"Cell volume {volume} is too small for rmax {rmax} (required >= {required_volume}).")

def compute_pairwise_rdf(atoms: Atoms,
                         ref_indices: List[int],
                         target_indices: List[int],
                         rmax: float, nbins: int,
                         volume: Optional[float] = None):
    """
    Compute pairwise radial distribution function between sets of atoms.
    
    Parameters
    ----------
    atoms : Atoms
        ASE Atoms object
    ref_indices : List[int]
        Indices of reference atoms
    target_indices : List[int]
        Indices of target atoms
    rmax : float
        Maximum radius for RDF calculation
    nbins : int
        Number of bins for histogram
    volume : float, optional
        System volume (calculated from atoms if None)
        
    Returns
    -------
    Tuple[np.ndarray, np.ndarray]
        RDF values and corresponding bin centers
    """
    N_total = len(atoms)
    dm = atoms.get_all_distances(mic=True)
    dr = float(rmax / nbins)
    volume = atoms.get_volume() if volume is None else volume

    
    if set(ref_indices) == set(target_indices):
        sub_dm = dm[np.ix_(ref_indices, target_indices)]
        sub_dm = np.triu(sub_dm, k=1)  # Exclude diagonal and use upper triangle
        distances = sub_dm[sub_dm > 0]
        
        N = len(ref_indices)
        
        # Division by 2 for same-species pairs
        norm = (4 * math.pi * dr * (N/volume) * N)/2
    else:
        sub_dm = dm[np.ix_(ref_indices, target_indices)]
        distances = sub_dm[sub_dm > 0]
        
        N_A = len(ref_indices)
        N_B = len(target_indices)
        
        norm = 4 * math.pi * dr * (N_A / volume) * N_total

    hist, bin_edges = np.histogram(distances, bins=nbins, range=(0, rmax))
    bin_centers = 0.5 * (bin_edges[1:] + bin_edges[:-1])
    rdf = hist / (norm * (bin_centers**2))
                  
    return rdf, bin_centers

class Analysis:
    """
    Class for analyzing atomic trajectories and calculating RDFs.
    
    Parameters
    ----------
    images : List[Atoms]
        List of ASE Atoms objects representing trajectory frames
    """
    
    def __init__(self, images: List[Atoms]):
        self.images = images

    def _get_slice(self, imageIdx: Optional[Union[int, slice]]):
        """
        Convert image index to slice for selecting trajectory frames.
        
        Parameters
        ----------
        imageIdx : Optional[Union[int, slice]]
            Index or slice to select images
            
        Returns
        -------
        slice
            Slice object for image selection
        """
        if imageIdx is None:
            return slice(None)
        return imageIdx

    def get_rdf(self, 
                rmax: float,
                nbins: int = 100,
                imageIdx: Optional[Union[int, slice]] = None,
                atomic_indices: Optional[Tuple[List[int], List[int]]] = None,
                return_dists: bool = False):
        """
        Calculate radial distribution function for trajectory frames.
        
        Parameters
        ----------
        rmax : float
            Maximum radius for RDF calculation
        nbins : int, optional
            Number of bins for histogram (default: 100)
        imageIdx : Optional[Union[int, slice]], optional
            Index or slice to select images (default: None, all images)
        atomic_indices : Optional[Tuple[List[int], List[int]]], optional
            Tuple of (reference_indices, target_indices) for partial RDF
        return_dists : bool, optional
            Whether to return bin center distances (default: False)
            
        Returns
        -------
        List[Union[np.ndarray, Tuple[np.ndarray, np.ndarray]]]
            List of RDF values or tuples (RDF, bin_centers) for each frame
        """
        sl = self._get_slice(imageIdx)
        images_to_process = self.images[sl]

        if atomic_indices is None:
            def process_image(image: Atoms):
                check_cell_and_r_max(image, rmax)
                full_indices = list(range(len(image)))
                rdf, bin_centers = compute_pairwise_rdf(image, full_indices, full_indices, rmax, nbins)
                return (rdf, bin_centers) if return_dists else rdf
        else:
            ref_indices, target_indices = atomic_indices
            def process_image(image: Atoms):
                check_cell_and_r_max(image, rmax)
                rdf, bin_centers = compute_pairwise_rdf(
                    image,  
                    ref_indices, 
                    target_indices,
                    rmax, 
                    nbins
                )
                return (rdf, bin_centers) if return_dists else rdf

        ls_rdf = Parallel(n_jobs=-1)(delayed(process_image)(image) for image in images_to_process)
        return ls_rdf

def analyze_rdf(use_prdf: bool,
                rmax: float,
                traj_path: str,
                nbins: int = 100,
                frame_skip: int = 10,
                output_filename: Optional[str] = None,
                atomic_indices: Optional[Tuple[List[int], List[int]]] = None,
                output_dir: str = 'custom_ase'):
    """
    Analyze trajectory and calculate radial distribution functions.
    
    Parameters
    ----------
    use_prdf : bool
        Whether to calculate partial RDF (True) or total RDF (False)
    rmax : float
        Maximum radius for RDF calculation
    traj_path : str
        Path to trajectory file
    nbins : int, optional
        Number of bins for histogram (default: 100)
    frame_skip : int, optional
        Number of frames to skip between analyses (default: 10)
    output_filename : Optional[str], optional
        Custom filename for output (default: None, auto-generated)
    atomic_indices : Optional[Tuple[List[int], List[int]]], optional
        Tuple of (reference_indices, target_indices) for partial RDF
    output_dir : str, optional
        Directory to save output files (default: 'custom_ase')
        
    Returns
    -------
    None
        Results are saved to file
        
    Raises
    ------
    ValueError
        If no images found in trajectory or if atomic_indices is missing for PRDF
    """
    images = read(traj_path, index=f'::{frame_skip}')
    if not images:
        raise ValueError("No images found in the trajectory.")
    # Check cell validity for the first image
    check_cell_and_r_max(images[0], rmax)
    
    analysis = Analysis(images)
    
    if use_prdf:
        if atomic_indices is None:
            raise ValueError("For partial RDF, atomic_indices must be provided.")
        ls_rdf = analysis.get_rdf(rmax, nbins, atomic_indices=atomic_indices, return_dists=True)
    else:
        ls_rdf = analysis.get_rdf(rmax, nbins, atomic_indices=None, return_dists=True)
    
    x_data_all = ls_rdf[0][1]
    y_data_all = [rdf for rdf, _ in ls_rdf]
    
    os.makedirs(output_dir, exist_ok=True)
    output_filename = output_filename or ('prdf_custom_indices.pkl' if use_prdf else 'rdf_total.pkl')
    output_filename = f'{output_filename}.pkl' if not output_filename.endswith('.pkl') else output_filename
    output_file = os.path.join(output_dir, output_filename)
    
    with open(output_file, 'wb') as f:
        pickle.dump({'x_data': x_data_all, 'y_data_all': y_data_all}, f)
    
    print(f"Data saved in '{output_file}'")
    return 

def main():
    """
    Parse command-line arguments and run RDF analysis.
    
    Command-line Parameters
    ----------------------
    --traj-path : str
        Path to trajectory file
    --rmax : float
        Maximum distance to consider for RDF
    --nbins : int
        Number of bins for the RDF histogram
    --frame-skip : int
        Use every n-th frame from trajectory
    --output-dir : str
        Directory to store output
    --output-filename : str
        Filename for output data
    --total-rdf
        Calculate total RDF (all atoms with all atoms)
    --partial-rdf
        Calculate partial RDF between specified atom groups
    --ref-indices : str
        Comma-separated list of reference atom indices for partial RDF
    --target-indices : str
        Comma-separated list of target atom indices for partial RDF
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('--traj-path', type=str, required=True)
    parser.add_argument('--rmax', type=float, required=True)
    parser.add_argument('--nbins', type=int, default=100)
    parser.add_argument('--frame-skip', type=int, default=10)
    parser.add_argument('--output-dir', type=str, default='custom_ase')
    parser.add_argument('--output-filename', type=str, default=None)
    
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--total-rdf', action='store_true')
    group.add_argument('--partial-rdf', action='store_true')
    
    parser.add_argument('--ref-indices', type=str, default=None)
    parser.add_argument('--target-indices', type=str, default=None)
    
    args = parser.parse_args()
    
    # Process atomic indices if partial RDF
    atomic_indices = None
    if args.partial_rdf:
        if not args.ref_indices or not args.target_indices:
            parser.error("--ref-indices and --target-indices must be provided with --partial-rdf")
        
        try:
            ref_indices = [int(i) for i in args.ref_indices.split(',')]
            target_indices = [int(i) for i in args.target_indices.split(',')]
            atomic_indices = (ref_indices, target_indices)
        except ValueError:
            parser.error("Indices must be comma-separated integers")
    
    # Run the analysis
    analyze_rdf(
        use_prdf=args.partial_rdf,
        rmax=args.rmax,
        traj_path=args.traj_path,
        nbins=args.nbins,
        frame_skip=args.frame_skip,
        output_filename=args.output_filename,
        atomic_indices=atomic_indices,
        output_dir=args.output_dir
    )

if __name__ == "__main__":
    main()

