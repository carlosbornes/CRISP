"""
Unit and regression test for the check_package package.
"""

# Import package, test suite, and other packages as needed
import sys

import pytest

import check_package


def test_check_package_imported():
    """Sample test, will always pass so long as import statement worked."""
    assert "check_package" in sys.modules
