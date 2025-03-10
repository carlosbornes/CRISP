import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import csv
from collections import defaultdict

def h_bond_2d(csv_file):
    # Read the CSV file and load the data into arrays
    data = pd.read_csv(csv_file)

    # Extract the columns from the data
    dist_oo = data['Dist(Donor-Acceptor)']  # Adjusted column name
    angle_ho1o2 = data['Angle(Hydrogen-Donor-Acceptor)']  # Adjusted column name

    # Set up the histogram bins
    dist_bins = np.linspace(dist_oo.min(), dist_oo.max(), 50)
    angle_bins = np.linspace(angle_ho1o2.min(), angle_ho1o2.max(), 50)

    # Create the 2D histogram
    hist, dist_edges, angle_edges = np.histogram2d(dist_oo, angle_ho1o2, bins=[dist_bins, angle_bins], weights=None)

    # Plot the 2D histogram
    fig, ax = plt.subplots()
    im = ax.imshow(hist.T, origin='lower', aspect='auto', extent=[dist_edges[0], dist_edges[-1], angle_edges[0], angle_edges[-1]])
    ax.set_xlabel('Distance (Donor-Acceptor)')
    ax.set_ylabel('Angle (H-Donor-Acceptor)')
    cbar = fig.colorbar(im)
    cbar.set_label('Counts')

    # Show the plot
    plt.show()

def plot_total_hydrogen_bonds_per_frame(csv_file):
    # Read the CSV file into a pandas DataFrame
    data = pd.read_csv(csv_file)

    # Extract frame indices and total hydrogen bond counts
    frame_indices = data['Frame']
    total_hydrogen_bonds = data['Total Hydrogen Bonds']

    # Create a line plot
    plt.figure(figsize=(10, 6))
    plt.plot(frame_indices, total_hydrogen_bonds, marker='o')
    plt.xlabel('Frame Index')
    plt.ylabel('Total Hydrogen Bonds')
    plt.title('Total Hydrogen Bonds per Frame')
    plt.grid(True)
    plt.show()

def count_double_donor_hydrogen_bonds(hydrogen_bonds_file):
    # Dictionary to store the count of unique frames
    frame_counts = defaultdict(int)
    
    # Read the hydrogen bonds data from the CSV file
    with open(hydrogen_bonds_file, 'r') as f:
        reader = csv.DictReader(f)
        donor_counts = defaultdict(int)
        current_frame = None
        
        for row in reader:
            frame = int(row['Frame'])
            donor = int(row['Donor'])  # Adjusted column name
            
            # Initialize the current frame if not set or if it has changed
            if frame != current_frame:
                if current_frame is not None:
                    # Calculate DD and SD counts for the current frame
                    dd_count = sum(1 for count in donor_counts.values() if count >= 2)
                    sd_count = sum(count for count in donor_counts.values()) - dd_count
                    frame_counts[current_frame] = {'DD': dd_count, 'SD': sd_count}
                
                current_frame = frame
                donor_counts.clear()
            
            # Count the occurrences of each donor atom for the current frame
            donor_counts[donor] += 1
    
    # Calculate counts for the last frame
    if current_frame is not None:
        dd_count = sum(1 for count in donor_counts.values() if count >= 2)
        sd_count = sum(count for count in donor_counts.values()) - dd_count
        frame_counts[current_frame] = {'DD': dd_count, 'SD': sd_count}
    
    return frame_counts

def write_results_to_file(output_file, frame_counts):
    # Write the results to the text file
    with open(output_file, 'w') as f:
        f.write("Frame\tSD\tDD\n")
        for frame, counts in frame_counts.items():
            dd_count = counts['DD']
            sd_count = counts['SD']
            f.write(f"{frame}\t{sd_count}\t{dd_count}\n")

    print(f"Results written to {output_file}")

def plot_hydrogen_bond_counts(text_file):
    # Read the data from the text file
    frame_numbers, sd_counts, dd_counts = [], [], []
    with open(text_file, 'r') as f:
        next(f)  # Skip the header line
        for line in f:
            frame, sd, dd = map(int, line.strip().split())
            frame_numbers.append(frame)
            sd_counts.append(sd)
            dd_counts.append(dd)
    
    # Create a plot to visualize the changes in SD and DD counts over time
    plt.figure(figsize=(10, 6))
    plt.plot(frame_numbers, sd_counts, label='SD Count', marker='o', linestyle='-', color='blue')
    plt.plot(frame_numbers, dd_counts, label='DD Count', marker='o', linestyle='-', color='red')
    plt.xlabel('Frame')
    plt.ylabel('Count')
    plt.title('Changes in SD and DD Hydrogen Bond Counts Over Time')
    plt.legend()
    plt.grid(True)
    
    # Show the plot
    plt.show()

if __name__ == "__main__":
    # Example usage in a script or Jupyter notebook

    # Perform hydrogen bond analysis and plot the results
    hydrogen_bonds_csv = 'hydrogen_bonds.csv'
    h_bond_2d(hydrogen_bonds_csv)

    # Plot total hydrogen bonds per frame
    total_hydrogen_bonds_csv = 'total_hydrogen_bonds_per_frame.csv'
    plot_total_hydrogen_bonds_per_frame(total_hydrogen_bonds_csv)
    
    hydrogen_bonds_csv = 'hydrogen_bonds.csv'
    # Count DD and SD hydrogen bonds for each frame
    frame_counts = count_double_donor_hydrogen_bonds(hydrogen_bonds_csv)
    hydrogen_bond_results_txt = 'hydrogen_bond_results.txt'

    # Write the results to the text file
    write_results_to_file(hydrogen_bond_results_txt, frame_counts)

    # Count double donor hydrogen bonds and plot the counts
    plot_hydrogen_bond_counts(hydrogen_bond_results_txt)
