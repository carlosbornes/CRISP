"""Command-line interface for CRISP package."""

import sys
import argparse
from CRISP._version import __version__
#from CRISP.tests.runner import run_tests

def main():
    """Main entry point for the CRISP command-line tool."""
    parser = argparse.ArgumentParser(description="CRISP: Concurrent Remote Interactive Simulation Program")
    parser.add_argument('--version', action='version', version=f'CRISP {__version__}')
        
    # Parse arguments
    parser.parse_args()
    
    if len(sys.argv) == 1:
        parser.print_help()
        
    return 0

if __name__ == "__main__":
    sys.exit(main())
