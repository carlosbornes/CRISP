# CRISP/data_analysis/contact.py

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import argparse
import os
from ase.io import read

def analyze_contacts(trajectory_file, type1_indices_file, type2_indices_file, 
                    cutoff_distance=3.0, tstep=50.0, skip=100, 
                    output_dir="contact_analysis", output_prefix="contacts_ow_ow",
                    plot_heatmap=True, plot_contact_counts=True):
    """
    Analyze contacts between two types of atoms in a molecular dynamics trajectory.
    
    Parameters:
        trajectory_file (str): Path to the ASE trajectory file
        type1_indices_file (str): Path to numpy file with type1 atom indices
        type2_indices_file (str): Path to numpy file with type2 atom indices
        cutoff_distance (float): Cutoff distance in Å to consider atoms in contact
        tstep (float): Time step per frame in fs
        skip (int): Process every 'skip' frame
        output_dir (str): Directory to save output files
        output_prefix (str): Prefix for output files
        plot_heatmap (bool): Whether to plot contact heatmap
        plot_contact_counts (bool): Whether to plot contact counts vs time
        
    Returns:
        dict: Dictionary containing analysis results
    """
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    print(f"Output files will be saved to {output_dir}")
    
    # Load atom indices
    type1_indices = np.load(type1_indices_file)
    type2_indices = np.load(type2_indices_file)
    
    print(f"Loaded {len(type1_indices)} type1 indices and {len(type2_indices)} type2 indices")
    
    # Keys are tuples (type1_index, type2_index)
    active_contacts = {}
    contact_times = {}

    # Initialize contact_times for each type1-type2 pair
    for type1 in type1_indices:
        for type2 in type2_indices:
            contact_times[(type1, type2)] = []

    # Define output file paths
    contacts_filename = os.path.join(output_dir, f"{output_prefix}.txt")
    stats_filename = os.path.join(output_dir, f"{output_prefix}_statistics.txt")
    matrix_csv_filename = os.path.join(output_dir, f"{output_prefix}_matrix.csv")
    
    # Lists for plotting contact counts vs time
    time_list = []
    contact_count_list = []
    
    # Process trajectory and calculate contacts
    with open(contacts_filename, "w") as contact_file:
        contact_file.write("Time_ps\t1_index\t2_index\tDistance_A\n")
        
        # Read the trajectory using ase.io.read
        traj = read(trajectory_file, ':')
        print(f"Processing trajectory with {len(traj)} frames, using every {skip}th frame")
        
        frame_index = 0
        for atoms in traj[::skip]:
            # Current simulation time in fs, then convert to ps:
            current_time_fs = frame_index * tstep * skip
            current_time_ps = current_time_fs / 1000.0
            time_list.append(current_time_ps)
            
            dist_matrix = atoms.get_all_distances(mic=True)        
            count_contacts = 0  # counter for the current frame
            
            # Loop over each type1 and type2 pair
            for type1 in type1_indices:
                for type2 in type2_indices:  # Fixed missing underscore
                    distance = dist_matrix[type1, type2]
                    if distance < cutoff_distance:
                        count_contacts += 1
                        # If this pair is entering the cutoff region, record the start time
                        if (type1, type2) not in active_contacts:
                            active_contacts[(type1, type2)] = current_time_ps
                        # Log this contact
                        contact_file.write(f"{current_time_ps:.3f}\t{type1}\t{type2}\t{distance:.3f}\n")
                    else:
                        # If the pair was active and now is no longer within cutoff
                        if (type1, type2) in active_contacts:
                            start_time = active_contacts.pop((type1, type2))
                            duration = current_time_ps - start_time
                            contact_times[(type1, type2)].append(duration)
            
            contact_count_list.append(count_contacts)
            frame_index += 1
            
            # Print progress every 10% of frames
            if frame_index % max(1, len(traj) // 10) == 0:
                print(f"Processed {frame_index} frames ({frame_index/len(traj)*100:.1f}%)")

    # After processing all frames, for any pairs still active, record the final event duration
    final_time_ps = frame_index * tstep * skip / 1000.0
    for pair, start_time in active_contacts.items():
        duration = final_time_ps - start_time
        contact_times[pair].append(duration)

    # Combine all contact event durations to compute the overall average contact time (in ps)
    all_event_times = []
    for events in contact_times.values():
        all_event_times.extend(events)

    results = {"time_list": time_list, "contact_count_list": contact_count_list}
    
    if all_event_times:
        overall_mean = np.mean(all_event_times)
        print(f"Overall average contact time: {overall_mean:.3f} ps")
        results["overall_mean_contact_time"] = overall_mean
    else:
        print("No contact events recorded.")
        results["overall_mean_contact_time"] = 0
    
    # Compute the mean contact time for a single type2 molecule
    if len(type2_indices) > 0:
        single_type2 = type2_indices[0]
        single_type2_events = []
        for type1 in type1_indices:
            single_type2_events.extend(contact_times.get((type1, single_type2), []))

        if single_type2_events:
            single_mean = np.mean(single_type2_events)
            print(f"Mean contact time for type2 {single_type2}: {single_mean:.3f} ps")
            results["single_type2_mean_contact_time"] = single_mean
        else:
            print(f"No contact events recorded for type2 {single_type2}.")
            results["single_type2_mean_contact_time"] = 0

    # Plot contact counts vs time if requested
    if plot_contact_counts:
        # Calculate the average contact count
        avg_contact_count = np.mean(contact_count_list)
        
        plt.figure(figsize=(10, 6))
        # Plot the contact counts
        plt.plot(time_list, contact_count_list, label='Contact counts')
        # Add a horizontal line for the average
        plt.axhline(y=avg_contact_count, color='r', linestyle='--', 
                   label=f'Average: {avg_contact_count:.2f}')
        
        plt.xlabel("Time (ps)")
        plt.ylabel("Number of Contacts")
        plt.title(f"Contact Counts Over Time (cutoff = {cutoff_distance} Å)")
        plt.grid(True)
        plt.legend()  # Add legend to show both lines
        plt.show()
        
        # Also save the plot
        counts_plot_filename = os.path.join(output_dir, f"{output_prefix}_counts_plot.png")
        plt.savefig(counts_plot_filename, dpi=300, bbox_inches='tight')
        print(f"Saved contacts count plot to {counts_plot_filename}")

    # Create and plot contact matrix if requested
    try:
        df = pd.read_csv(contacts_filename, sep="\t")
        
        # Get unique indices
        all_1_indices = pd.Index(df['1_index'].unique())
        all_2_indices = pd.Index(df['2_index'].unique())
        
        # Reindex the contact matrix
        contact_matrix = df.groupby(['1_index', '2_index']).size().unstack(fill_value=0)
        contact_matrix = contact_matrix.reindex(index=all_1_indices, columns=all_2_indices, fill_value=0)
        
        # Save contact matrix to CSV
        contact_matrix.to_csv(matrix_csv_filename)
        print(f"Saved contact matrix to {matrix_csv_filename}")
        
        results["contact_matrix"] = contact_matrix
        
        if plot_heatmap:
            # Plot the contact matrix
            plt.figure(figsize=(14, 12))
            heatmap = sns.heatmap(
                contact_matrix, 
                annot=False,
                cmap="viridis", 
                cbar_kws={'label': 'Number of Contacts'},
                xticklabels=contact_matrix.columns, 
                yticklabels=contact_matrix.index
            )
            
            plt.xlabel("Type 2 Indices", fontsize=12)
            plt.ylabel("Type 1 Indices", fontsize=12)
            plt.title(f"Contact Matrix (cutoff = {cutoff_distance} Å)", fontsize=16)
            plt.xticks(rotation=90)
            plt.yticks(rotation=0)
            plt.tight_layout()
            
            # Save the heatmap before showing it
            heatmap_filename = os.path.join(output_dir, f"{output_prefix}_heatmap.png")
            plt.savefig(heatmap_filename, dpi=300, bbox_inches='tight')
            print(f"Saved heatmap to {heatmap_filename}")
            
            plt.show()
        
        # Calculate statistics for writing to file
        flattened_contacts = contact_matrix.stack()
        
        # Top 10 contact pairs
        top_10_contacts = flattened_contacts.sort_values(ascending=False).head(10)
        
        # Bottom 10 non-zero contact pairs
        nonzero_contacts = flattened_contacts[flattened_contacts > 0]
        bottom_10_contacts = nonzero_contacts.sort_values().head(10)
        
        # Calculate total contacts per type
        total_contacts_per_type1 = contact_matrix.sum(axis=1).sort_values(ascending=False)
        total_contacts_per_type2 = contact_matrix.sum(axis=0).sort_values(ascending=False)
        
        # Calculate lowest non-zero contacts per type
        nonzero_type1 = total_contacts_per_type1[total_contacts_per_type1 > 0]
        nonzero_type2 = total_contacts_per_type2[total_contacts_per_type2 > 0]
        bottom_type1 = nonzero_type1.sort_values().head(10)
        bottom_type2 = nonzero_type2.sort_values().head(10)
        
        # Write statistics to file
        with open(stats_filename, "w") as stats_file:
            stats_file.write("Contact Analysis Statistics\n")
            stats_file.write(f"=========================\n\n")
            
            # Top 10 Type1-Type2 Contact Pairs
            stats_file.write("Top 10 Type1-Type2 Contact Pairs (Highest Counts):\n")
            stats_file.write("------------------------------------------------\n")
            for (i, j), count in top_10_contacts.items():
                stats_file.write(f"Indices ({i}, {j}) - Contacts: {count}\n")
            stats_file.write("\n")
            
            # Bottom 10 Type1-Type2 Contact Pairs (non-zero)
            stats_file.write("Bottom 10 Type1-Type2 Contact Pairs (Non-Zero):\n")
            stats_file.write("-----------------------------------------------\n")
            for (i, j), count in bottom_10_contacts.items():
                stats_file.write(f"Indices ({i}, {j}) - Contacts: {count}\n")
            stats_file.write("\n")
            
            # Top 10 Type1 Indices by Total Contacts
            stats_file.write("Top 10 Type1 Indices by Total Contacts:\n")
            stats_file.write("--------------------------------------\n")
            for index, total in total_contacts_per_type1.head(10).items():
                stats_file.write(f"Index {index} - Total Contacts: {total}\n")
            stats_file.write("\n")
            
            # Top 10 Type2 Indices by Total Contacts
            stats_file.write("Top 10 Type2 Indices by Total Contacts:\n")
            stats_file.write("--------------------------------------\n")
            for index, total in total_contacts_per_type2.head(10).items():
                stats_file.write(f"Index {index} - Total Contacts: {total}\n")
            stats_file.write("\n")
            
            # Bottom 10 Type1 Indices by Total Contacts (non-zero)
            stats_file.write("Bottom 10 Type1 Indices by Total Contacts (Non-Zero):\n")
            stats_file.write("--------------------------------------------------\n")
            for index, total in bottom_type1.items():
                stats_file.write(f"Index {index} - Total Contacts: {total}\n")
            stats_file.write("\n")
            
            # Bottom 10 Type2 Indices by Total Contacts (non-zero)
            stats_file.write("Bottom 10 Type2 Indices by Total Contacts (Non-Zero):\n")
            stats_file.write("--------------------------------------------------\n")
            for index, total in bottom_type2.items():
                stats_file.write(f"Index {index} - Total Contacts: {total}\n")
        
        print(f"Saved contact statistics to {stats_filename}")
            
        results["top_10_contacts"] = top_10_contacts
        results["bottom_10_contacts"] = bottom_10_contacts
        results["total_contacts_per_type1"] = total_contacts_per_type1
        results["total_contacts_per_type2"] = total_contacts_per_type2
        results["bottom_type1_contacts"] = bottom_type1
        results["bottom_type2_contacts"] = bottom_type2
        
    except Exception as e:
        print(f"Error analyzing contact data: {e}")
    
    return results

def main():
    """Main function to parse arguments and run contact analysis."""
    parser = argparse.ArgumentParser(
        description="Analyze contacts between two types of atoms in a molecular dynamics trajectory."
    )
    
    parser.add_argument(
        "trajectory",
        type=str,
        help="Path to the ASE trajectory file."
    )
    
    parser.add_argument(
        "--type1_indices",
        type=str,
        required=True,
        help="Path to numpy file with type1 atom indices."
    )
    
    parser.add_argument(
        "--type2_indices",
        type=str,
        required=True,
        help="Path to numpy file with type2 atom indices."
    )
    
    parser.add_argument(
        "--cutoff",
        type=float,
        default=3.0,
        help="Cutoff distance in Å to consider atoms in contact. Default is 3.0 Å."
    )
    
    parser.add_argument(
        "--tstep",
        type=float,
        default=50.0,
        help="Time step per frame in fs. Default is 50.0 fs."
    )
    
    parser.add_argument(
        "--skip",
        type=int,
        default=100,
        help="Process every 'skip' frame. Default is 100."
    )
    
    parser.add_argument(
        "--output_dir",
        type=str,
        default="contact_analysis",
        help="Directory to save output files. Default is 'contact_analysis'."
    )
    
    parser.add_argument(
        "--output_prefix",
        type=str,
        default="contacts_ow_ow",
        help="Prefix for output files. Default is 'contacts_ow_ow'."
    )
    
    parser.add_argument(
        "--no_heatmap",
        action="store_true",
        help="Disable heatmap plot display."
    )
    
    parser.add_argument(
        "--no_counts_plot",
        action="store_true",
        help="Disable contact counts plot display."
    )
    
    args = parser.parse_args()
    
    # Run contact analysis
    analyze_contacts(
        trajectory_file=args.trajectory,
        type1_indices_file=args.type1_indices,
        type2_indices_file=args.type2_indices,
        cutoff_distance=args.cutoff,
        tstep=args.tstep,
        skip=args.skip,
        output_dir=args.output_dir,
        output_prefix=args.output_prefix,
        plot_heatmap=not args.no_heatmap,
        plot_contact_counts=not args.no_counts_plot
    )

if __name__ == "__main__":
    main()

