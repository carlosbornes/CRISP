# CRISP/tests/test_CRISP.py

import unittest
import os
from CRISP.atom_indices import TrajectoryAnalyzer

class TestTrajectoryAnalyzer(unittest.TestCase):
    def setUp(self):
        # Use a path to a test trajectory file that you know exists
        self.test_traj_path = 'path/to/test_traj.traj'
        self.analyzer = TrajectoryAnalyzer(self.test_traj_path)

    def test_find_atom_indices(self):
        indices = self.analyzer.find_atom_indices()
        self.assertEqual(len(indices), 7)  # Adjust this based on your expected output

    def test_save_indices(self):
        self.analyzer.save_indices()
        output_folder = os.path.join(os.path.dirname(self.test_traj_path), 'indices')
        self.assertTrue(os.path.exists(os.path.join(output_folder, 'lengths.npy')))

if __name__ == '__main__':
    unittest.main()
