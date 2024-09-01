from ase.io import Trajectory
import numpy as np
from joblib import Parallel, delayed
import pickle
import argparse

# Define a function to calculate the distance between two atoms
def distance(atoms, i, j):
    return atoms.get_distance(i, j, mic=True)

def classify_atom_coordination(structure, atom_indices, symbol_pairs, cutoffs):
    coordination_types = {}

    for atom_index in atom_indices:
        coordination = 0
        atom_symbol = structure[atom_index].symbol
        for other_index, other_atom in enumerate(structure):
            if other_index != atom_index:
                other_symbol = other_atom.symbol
                if (atom_symbol, other_symbol) in cutoffs:
                    dist = distance(structure, atom_index, other_index)
                    if dist < cutoffs[(atom_symbol, other_symbol)]:
                        coordination += 1
                elif (other_symbol, atom_symbol) in cutoffs:
                    dist = distance(structure, atom_index, other_index)
                    if dist < cutoffs[(other_symbol, atom_symbol)]:
                        coordination += 1

        coordination_types[atom_index] = coordination

    return coordination_types

def analyze_atom_coordination(atoms, atom_indices, symbol_pairs, cutoffs):
    atom_coordination_types = classify_atom_coordination(atoms, atom_indices, symbol_pairs, cutoffs)
    return atom_coordination_types

def calculate_atom_coordination(file_path, indices_path, output_file, frame_skip, cutoffs, atom):
    # Load the trajectory
    traj = Trajectory(file_path)
    frames = [frame for i, frame in enumerate(traj) if i % frame_skip == 0]

    # Load custom indices
    custom_indices = np.load(indices_path)

    # Filter the indices based on the atom symbol
    custom_indices = [i for i in custom_indices if traj[0][i].symbol == atom]

    # Use Parallel from joblib to analyze frames in parallel
    coordination_types_list = Parallel(n_jobs=-1)(
        delayed(analyze_atom_coordination)(atoms, custom_indices, cutoffs.keys(), cutoffs) for atoms in frames
    )

    # Save the coordination data for further analysis
    with open(output_file, 'wb') as f:
        pickle.dump(coordination_types_list, f)

    # For testing purposes, print the first entry of the coordination types list
    print("First entry of the coordination types list:", coordination_types_list[0])


def parse_and_run():
    parser = argparse.ArgumentParser(description="Calculate atom coordination in a trajectory.")
    parser.add_argument('trajectory', type=str, help='Path to ASE trajectory file (.traj)')
    parser.add_argument('indices', type=str, help='Path to numpy array containing custom indices')
    parser.add_argument('output', type=str, default='coordination_data.pkl', help='Output file name for coordination data')
    parser.add_argument('--frame_skip', type=int, default=10, help='Number of frames to skip. Default is 10.')
    parser.add_argument('--cutoffs', type=str, required=True, help='Cutoff distances for atom pairs in the format symbol1-symbol2:cutoff,symbol3-symbol4:cutoff')
    parser.add_argument('--atom', type=str, required=True, help='Symbol of the atom to check the coordination for.')
    
    args = parser.parse_args()

    cutoffs = {}
    symbol_pairs = []
    pairs = args.cutoffs.split(',')
    for pair in pairs:
        symbols, cutoff = pair.split(':')
        symbol1, symbol2 = symbols.split('-')
        cutoffs[(symbol1, symbol2)] = float(cutoff)
        symbol_pairs.append((symbol1, symbol2))

    calculate_atom_coordination(args.trajectory, args.indices, args.output, args.frame_skip, symbol_pairs, cutoffs)

if __name__ == "__main__":
    parse_and_run()
