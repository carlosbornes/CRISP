#!/usr/bin/env python

"""
Test runner for the CRISP package.
Run this with: python run_tests.py
"""

import os
import sys
import unittest
import subprocess

def run_tests():
    """Run the CRISP test suite."""
    print("Running CRISP test suite...")
    
    # Get the project root directory (where this script is located)
    project_root = os.path.dirname(os.path.abspath(__file__))
    
    # Get the tests directory
    test_dir = os.path.join(project_root, "CRISP", "tests")
    
    # First run unittest-based tests
    print("\n=== Running unittest tests ===")
    test_loader = unittest.TestLoader()
    test_suite = test_loader.discover(test_dir, pattern="test_*.py")
    test_runner = unittest.TextTestRunner(verbosity=2)
    unittest_result = test_runner.run(test_suite)
    
    # Then run pytest-based tests if pytest is available
    print("\n=== Running pytest tests ===")
    try:
        import pytest
        pytest_args = [test_dir, "-v"]
        pytest_result = pytest.main(pytest_args)
        if pytest_result != 0:
            print("Pytest tests failed")
            return 1
    except ImportError:
        print("Pytest not available. Skipping pytest tests.")
    
    if not unittest_result.wasSuccessful():
        print("Unittest tests failed")
        return 1
    
    print("\nAll tests completed successfully!")
    return 0

if __name__ == "__main__":
    sys.exit(run_tests())
