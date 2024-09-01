import os
import numpy as np
import ase.io
import csv


def atom_indices(file_path, frame_index=0, cutoffs=None):
    new = ase.io.read(file_path, index=frame_index, format="traj")
    dist_matrix = new.get_all_distances(mic=True)
    symbols = new.get_chemical_symbols()

    unique_symbols = list(set(symbols))

    indices_by_symbol = {symbol: [] for symbol in unique_symbols}

    for idx, atom in enumerate(new):
        indices_by_symbol[atom.symbol].append(idx)

    cutoff_indices = {}

    if cutoffs:
        for pair, cutoff in cutoffs.items():
            symbol1, symbol2 = pair
            pair_indices_distances = []
            if symbol1 in indices_by_symbol and symbol2 in indices_by_symbol:
                for idx1 in indices_by_symbol[symbol1]:
                    for idx2 in indices_by_symbol[symbol2]:
                        if dist_matrix[idx1, idx2] < cutoff:
                            pair_indices_distances.append(
                                (idx1, idx2, dist_matrix[idx1, idx2])
                            )
            cutoff_indices[pair] = pair_indices_distances

    return indices_by_symbol, dist_matrix, cutoff_indices


def run_atom_indices(file_path, output_folder, frame_index=0, cutoffs=None):
    indices, dist_matrix, cutoff_indices = atom_indices(file_path, frame_index, cutoffs)

    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    lengths = {symbol: len(indices[symbol]) for symbol in indices}
    np.save(os.path.join(output_folder, "lengths.npy"), lengths)

    for symbol, data in indices.items():
        np.save(os.path.join(output_folder, f"{symbol}_indices.npy"), data)
        print(f"Length of {symbol} indices: {len(data)}")

    print("Outputs saved.")

    # Save cutoff indices into CSV files
    cutoff_folder = os.path.join(output_folder, "cutoff")
    if not os.path.exists(cutoff_folder):
        os.makedirs(cutoff_folder)

    for pair, pair_indices_distances in cutoff_indices.items():
        symbol1, symbol2 = pair
        filename = f"{symbol1}-{symbol2}_cutoff.csv"
        filepath = os.path.join(cutoff_folder, filename)
        with open(filepath, mode="w", newline="") as file:
            writer = csv.writer(file)
            writer.writerow([f"{symbol1} index", f"{symbol2} index", "distance"])
            writer.writerows(pair_indices_distances)
        print(f"Saved cutoff indices for {symbol1}-{symbol2} to {filepath}")


def main():
    import argparse

    parser = argparse.ArgumentParser(
        description="Analyze ASE trajectory and extract atom indices by symbol."
    )
    parser.add_argument("file_path", type=str, help="Path to the ASE trajectory file.")
    parser.add_argument(
        "output_folder", type=str, help="Folder to save the output indices."
    )
    parser.add_argument(
        "--frame_index",
        type=int,
        default=0,
        help="Frame index of the trajectory file to analyze. Default is 0.",
    )
    parser.add_argument(
        "--cutoffs",
        type=str,
        help="Cutoff distances for atomic symbol pairs in the format symbol1-symbol2:cutoff,symbol3-symbol4:cutoff",
        default=None,
    )

    args = parser.parse_args()

    cutoffs = None
    if args.cutoffs:
        cutoffs = {}
        pairs = args.cutoffs.split(",")
        for pair in pairs:
            symbols, cutoff = pair.split(":")
            symbol1, symbol2 = symbols.split("-")
            cutoffs[(symbol1, symbol2)] = float(cutoff)

    run_atom_indices(args.file_path, args.output_folder, args.frame_index, cutoffs)


if __name__ == "__main__":
    main()
