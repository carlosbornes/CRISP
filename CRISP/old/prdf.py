import os
import numpy as np
import pickle
from ase.io import read
from ase import Atoms
from typing import Optional, Union
from joblib import Parallel, delayed

def compute_pairwise_rdf(atoms, ref_indices, target_indices, rmax, nbins, volume=None):
    distances = []
    cell = atoms.get_cell()
    for ref_idx in ref_indices:
        for target_idx in target_indices:
            if ref_idx != target_idx:
                dist = atoms.get_distance(ref_idx, target_idx, mic=True)
                distances.append(dist)
    
    distances = np.array(distances)
    hist, bin_edges = np.histogram(distances, bins=nbins, range=(0, rmax))
    bin_centers = 0.5 * (bin_edges[1:] + bin_edges[:-1])
    rdf = hist / (4 * np.pi * bin_centers**2 * (rmax / nbins))
    
    if volume:
        rdf /= volume
    else:
        vol = atoms.get_volume()
        rdf /= vol
    
    return rdf, bin_centers

class Analysis:
    def __init__(self, images):
        self.images = images

    def _get_slice(self, imageIdx):
        if imageIdx is None:
            return slice(None)
        return imageIdx

    def get_rdf(self, rmax: float, nbins: int = 100, imageIdx: Union[int, slice, None] = None,
                elements: Optional[list] = None, return_dists: bool = False, volume: Optional[float] = None,
                pairwise: Optional[tuple] = None):
        """Get RDF.

        Parameters:
        rmax: float
            Maximum distance of RDF.
        nbins: int
            Number of bins to divide RDF.
        imageIdx: int/slice/None
            Images to analyze, see :func:`_get_slice` for details.
        elements: list
            List of indices to consider for RDF calculation.
        pairwise: tuple
            Calculate RDF for specific pair interactions (e.g., ('Al', 'O')).

        Returns:
        return: list of lists / list of tuples of lists
            If return_dists is True, the returned tuples contain (rdf, distances). Otherwise
            only rdfs for each image are returned.
        """

        sl = self._get_slice(imageIdx)
        images_to_process = self.images[sl]

        if elements is None:
            raise ValueError("Elements indices must be specified for RDF calculation!")
        if pairwise is None:
            raise ValueError("Pairwise interactions must be specified for RDF calculation!")
        
        def process_image(image):
            tmp_image = Atoms(cell=image.get_cell(), pbc=image.get_pbc())
            original_to_tmp_idx = {}
            for idx in elements:
                original_to_tmp_idx[idx] = len(tmp_image)
                tmp_image.append(image[idx])

            ref_idxs = [original_to_tmp_idx[idx] for idx in elements if image[idx].symbol == pairwise[0]]
            target_idxs = [original_to_tmp_idx[idx] for idx in elements if image[idx].symbol == pairwise[1]]

            rdf, distances = compute_pairwise_rdf(tmp_image, ref_idxs, target_idxs, rmax, nbins, volume)
            return (rdf, distances) if return_dists else rdf
        
        ls_rdf = Parallel(n_jobs=-1)(delayed(process_image)(image) for image in images_to_process)

        return ls_rdf

def analyze_rdf(pairwise, rmax, traj_path, atom1_array_path, atom2_array_path, nbins=100, use_prdf=True, frame_skip=10, output_filename=None):
    # Load the trajectory
    images = read(traj_path, index='::{}'.format(frame_skip))

    if use_prdf:
        # Load custom indices from NumPy arrays
        atom1_array = np.load(atom1_array_path)
        if pairwise[0] == pairwise[1]:
            combined_array = atom1_array
        else:
            atom2_array = np.load(atom2_array_path)
            combined_array = np.concatenate([atom1_array, atom2_array])
    else:
        # Extract indices based on atomic symbols
        combined_array = []
        for image in images:
            for idx, atom in enumerate(image):
                if atom.symbol in pairwise:
                    combined_array.append(idx)
        combined_array = np.array(combined_array)

    # Initialize the Analysis class with the loaded images
    analysis = Analysis(images)

    # Get RDF for the provided indices
    rdf_results = analysis.get_rdf(rmax, nbins, elements=combined_array.tolist(), return_dists=True, volume=True, pairwise=tuple(pairwise))

    # Extract RDF data and distances
    x_data_all = np.linspace(0, rmax, nbins)
    y_data_all = []
    x_data = []
    y_data = []

    for rdf, distances in rdf_results:
        x_data.append(distances)
        y_data.append(rdf)
        y_data_all.append(rdf)

    # Ensure the PRDF directory exists
    output_dir = 'PRDF'
    os.makedirs(output_dir, exist_ok=True)

    # Determine the output file name
    if output_filename is None:
        output_filename = f'rdf_{"_".join(pairwise)}.pkl'
    else:
        output_filename = f'{output_filename}.pkl'
    
    output_file = os.path.join(output_dir, output_filename)

    # Save the data in a pickle file
    data = {'x_data': x_data, 'y_data': y_data, 'x_data_all': x_data_all, 'y_data_all': y_data_all}
    with open(output_file, 'wb') as f:
        pickle.dump(data, f)

    # Print to verify
    print(f"Data saved in '{output_file}'")

    return x_data_all, y_data_all


