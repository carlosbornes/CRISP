from ase.io.trajectory import Trajectory
import numpy as np
import csv
from joblib import Parallel, delayed
import argparse

# Define a function to calculate the distance between two atoms
def distance(atoms, i, j):
    return atoms.get_distance(i, j, mic=True)

def count_hydrogen_bonds(atoms, frame_idx, donor_indices, acceptor_indices, hydrogen_indices, donor_acceptor_distance_threshold, donor_hydrogen_distance_threshold, angle_cutoff):
    # Loop over all pairs of donor and acceptor atoms and check whether there is a hydrogen atom between them
    # Create an empty list to store the data
    data = []
    num_hydrogen_bonds = 0  # Initialize the count for total hydrogen bonds
    
    for donor in donor_indices:
        for acceptor in acceptor_indices:
            # Check whether there is a hydrogen atom between the donor and acceptor atoms
            h_between = False
            for hydrogen in hydrogen_indices:
                if distance(atoms, hydrogen, donor) < donor_hydrogen_distance_threshold and distance(atoms, donor, acceptor) < donor_acceptor_distance_threshold:
                    try:
                        angle = atoms.get_angle(hydrogen, donor, acceptor, mic=True)
                        if angle < angle_cutoff:
                            h_between = True
                            break
                    except ZeroDivisionError:
                        # Skip this pair if angle calculation results in ZeroDivisionError
                        continue
            # If there is a hydrogen atom between the donor and acceptor atoms and the angle is below the cutoff,
            # append the data to the list and increment the count
            if h_between:
                data.append([frame_idx, donor, hydrogen, acceptor, distance(atoms, hydrogen, donor), distance(atoms, donor, acceptor), angle])
                num_hydrogen_bonds += 1

    return data, num_hydrogen_bonds

def analyze_hydrogen_bonds(traj_path, donor_indices_path, acceptor_indices_path, hydrogen_indices_path, output_filename='hydrogen_bonds.csv', frame_skip=1, donor_acceptor_distance_threshold=3.5, donor_hydrogen_distance_threshold=1.2, angle_cutoff=30.0):
    traj = Trajectory(traj_path)
    donor_indices = np.load(donor_indices_path)
    acceptor_indices = np.load(acceptor_indices_path)
    hydrogen_indices = np.load(hydrogen_indices_path)

    # Parallelize the hydrogen bond analysis using joblib
    all_data = Parallel(n_jobs=-1)(delayed(count_hydrogen_bonds)(
        atoms, frame_idx, donor_indices, acceptor_indices, hydrogen_indices,
        donor_acceptor_distance_threshold, donor_hydrogen_distance_threshold, angle_cutoff
    ) for frame_idx, atoms in enumerate(traj[::frame_skip]))

    data = [entry for sublist, _ in all_data for entry in sublist]
    total_hydrogen_bonds_per_frame = [num_bonds for _, num_bonds in all_data]

    with open(output_filename, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['Frame', 'Donor', 'Hydrogen', 'Acceptor', 'Dist(Donor-Hydrogen)', 'Dist(Donor-Acceptor)', 'Angle(Hydrogen-Donor-Acceptor)'])
        writer.writerows(data)

    with open('total_hydrogen_bonds_per_frame.csv', 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['Frame', 'Total Hydrogen Bonds'])
        for frame_idx, num_bonds in enumerate(total_hydrogen_bonds_per_frame):
            writer.writerow([frame_idx, num_bonds])

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Hydrogen bond analysis.")
    parser.add_argument('--trajectory', type=str, help='Path to ASE trajectory file (.traj)')
    parser.add_argument('--donor_indices', type=str, help='Path to numpy array containing donor indices')
    parser.add_argument('--acceptor_indices', type=str, help='Path to numpy array containing acceptor indices')
    parser.add_argument('--hydrogen_indices', type=str, help='Path to numpy array containing hydrogen indices')
    parser.add_argument('--output', type=str, default='hydrogen_bonds.csv', help='Output file name for hydrogen bond data')
    parser.add_argument('--frame_skip', type=int, default=1, help='Number of frames to skip while processing the trajectory')
    parser.add_argument('--donor_acceptor_distance_threshold', type=float, default=3.5, help='Donor-acceptor distance threshold (default: 3.5 Å)')
    parser.add_argument('--donor_hydrogen_distance_threshold', type=float, default=1.2, help='Donor-hydrogen distance threshold (default: 1.2 Å)')
    parser.add_argument('--angle_cutoff', type=float, default=30.0, help='Angle cutoff for hydrogen bonds (default: 30.0 degrees)')
    args = parser.parse_args()
    
    analyze_hydrogen_bonds(
        args.trajectory, 
        args.donor_indices, 
        args.acceptor_indices, 
        args.hydrogen_indices,
        args.output, 
        args.frame_skip, 
        args.donor_acceptor_distance_threshold, 
        args.donor_hydrogen_distance_threshold, 
        args.angle_cutoff
    )
