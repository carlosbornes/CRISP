import numpy as np
import ase.io
from ase.visualize import view


def visualize(traj_file, frame_index=0):
    # Load atoms from .traj file
    atoms = ase.io.read(traj_file, index=frame_index, format="traj")
    
    view(atoms)


def main():
    import argparse
    parser = argparse.ArgumentParser(description='Visualize atom indices from a .traj file.')
    parser.add_argument('traj_file', type=str, help='Path to the ASE trajectory (.traj) file.')
#    parser.add_argument('indices_list', type=int, nargs='+', help='List of atom indices to visualize.')
    parser.add_argument('--frame_index', type=int, default=0, help='Frame index of the trajectory file to analyze. Default is 0.')
    args = parser.parse_args()

    # Visualize atom indices
    visualize(args.trag_file, args.frame_index)

if __name__ == "__main__":
    main()
