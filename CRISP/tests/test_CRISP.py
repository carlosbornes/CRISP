import unittest
import os
import numpy as np
from CRISP.atom_indices import atom_indices, run_atom_indices
import shutil

class TestAtomIndices(unittest.TestCase):
    def setUp(self):
        # Update the path to reflect the correct location of the test file
        self.test_traj_path = '/mnt/c/Users/sahaC/Desktop/side_project/molssi_best_practices/CRISP/CRISP/data/wrapped_traj.traj'
        self.output_folder = '/mnt/c/Users/sahaC/Desktop/side_project/molssi_best_practices/CRISP/CRISP/data/indices'
        
        if not os.path.exists(self.output_folder):
            os.makedirs(self.output_folder)

    def test_atom_indices(self):
        indices, dist_matrix, cutoff_indices = atom_indices(self.test_traj_path)
        
        # Example checks - adjust based on your expectations
        self.assertIsInstance(indices, dict)
        self.assertIsInstance(dist_matrix, np.ndarray)
        self.assertIsInstance(cutoff_indices, dict)
        self.assertIn('H', indices)  # Adjust based on your specific data

    def test_run_atom_indices(self):
        run_atom_indices(self.test_traj_path, self.output_folder)
        
        # Check if output files are created
        lengths_path = os.path.join(self.output_folder, 'lengths.npy')
        self.assertTrue(os.path.exists(lengths_path))

        # Check if length of indices is correct
        lengths = np.load(lengths_path, allow_pickle=True).item()
        self.assertIn('H', lengths)  # Adjust based on your specific data
        self.assertGreater(lengths['H'], 0)  # Adjust based on your specific data

def tearDown(self):
    # Clean up output files and directories
    if os.path.exists(self.output_folder):
        for file in os.listdir(self.output_folder):
            file_path = os.path.join(self.output_folder, file)
            if os.path.isfile(file_path):
                os.remove(file_path)
            elif os.path.isdir(file_path):
                shutil.rmtree(file_path)  # Use shutil to remove directories
        os.rmdir(self.output_folder)


if __name__ == '__main__':
    unittest.main()
