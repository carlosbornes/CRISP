import unittest
import os
import glob
import numpy as np
from CRISP.simulation_utility.atomic_indices import atom_indices, run_atom_indices
import shutil

class TestAtomIndices(unittest.TestCase):
    def setUp(self):
        # Get the directory of the test file
        test_dir = os.path.dirname(__file__)
        
        # Navigate one level up to locate the 'data' directory
        data_dir = os.path.abspath(os.path.join(test_dir, '..', 'data'))
        
        # Find all .traj files in the data directory
        traj_files = glob.glob(os.path.join(data_dir, '*.traj'))
        
        if not traj_files:
            self.fail(f"No .traj file found in {data_dir}.")
        
        # Use the first .traj file found
        self.test_traj_path = traj_files[0]
        self.output_folder = os.path.join(data_dir, 'indices')

        if not os.path.exists(self.output_folder):
            os.makedirs(self.output_folder)

    def test_atom_indices(self):
        indices, dist_matrix, cutoff_indices = atom_indices(self.test_traj_path)
        
        self.assertIsInstance(indices, dict)
        self.assertIsInstance(dist_matrix, np.ndarray)
        self.assertIsInstance(cutoff_indices, dict)
        self.assertIn('H', indices)

    def test_run_atom_indices(self):
        run_atom_indices(self.test_traj_path, self.output_folder)
        
        lengths_path = os.path.join(self.output_folder, 'lengths.npy')
        self.assertTrue(os.path.exists(lengths_path))

        lengths = np.load(lengths_path, allow_pickle=True).item()
        self.assertIn('H', lengths)
        self.assertGreater(lengths['H'], 0)

    def tearDown(self):
        # Clean up output files and directories
        if os.path.exists(self.output_folder):
            for file in os.listdir(self.output_folder):
                file_path = os.path.join(self.output_folder, file)
                if os.path.isfile(file_path):
                    os.remove(file_path)
                elif os.path.isdir(file_path):
                    shutil.rmtree(file_path)
            os.rmdir(self.output_folder)

if __name__ == '__main__':
    unittest.main()
