# CRISP/oxygen_coordination_analysis.py

import pickle
import matplotlib.pyplot as plt
import numpy as np
import itertools

def load_coordination_data(pickle_file):
    with open(pickle_file, 'rb') as f:
        return pickle.load(f)

def calculate_avg_percentages(coordination_data):
    # Get unique coordination types from the data
    coord_types = set(itertools.chain.from_iterable(frame_data.values() for frame_data in coordination_data))
    coord_types = sorted(coord_types)  # Sort the coordination types for consistent plotting
    avg_percentages = {coord_type: [] for coord_type in coord_types}

    for frame_data in coordination_data:
        total_oxygen_atoms = len(frame_data)
        
        for coord_type in coord_types:
            count = sum(1 for value in frame_data.values() if value == coord_type)
            avg_percentage = count / total_oxygen_atoms * 100
            avg_percentages[coord_type].append(avg_percentage)
    
    return avg_percentages

def plot_avg_percentages(avg_percentages, output_file=None):
    frames = list(range(len(next(iter(avg_percentages.values())))))
    coord_types = list(avg_percentages.keys())
    
    # Generate colors and markers based on the number of coordination types
    colors = plt.cm.get_cmap('tab10', len(coord_types)).colors
    markers = itertools.cycle(['o', 's', 'D', '^', 'v', 'p', '*', '+', 'x'])

    plt.figure(figsize=(10, 6))

    for i, coord_type in enumerate(coord_types):
        plt.plot(frames, avg_percentages[coord_type], label=f'{coord_type}-coordinated', color=colors[i], marker=next(markers))

    plt.xlabel('Frame')
    plt.ylabel('Average Percentage')
    plt.title('Average Percentage of Coordinated Atoms')
    
    # Place legend outside the plot
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    
    plt.grid(True)
    
    if output_file:
        plt.savefig(output_file, bbox_inches='tight')  # Use bbox_inches='tight' to include legend outside
    else:
        plt.show()

def calculate_overall_avg_percentages(avg_percentages):
    return {coord_type: sum(percentages) / len(percentages) for coord_type, percentages in avg_percentages.items()}

def main(pickle_file, plot_output_file=None):
    coordination_data = load_coordination_data(pickle_file)
    avg_percentages = calculate_avg_percentages(coordination_data)
    plot_avg_percentages(avg_percentages, plot_output_file)
    overall_avg_percentages = calculate_overall_avg_percentages(avg_percentages)
    
    # Print the overall average percentages
    for coord_type, avg_percentage in overall_avg_percentages.items():
        print(f'Overall Average percentage of {coord_type}-coordinated atoms: {avg_percentage:.2f}%')

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Analyze and plot coordination data.")
    parser.add_argument('pickle_file', type=str, help='Path to the pickle file containing coordination data')
    parser.add_argument('--plot_output', type=str, default=None, help='Path to save the plot (optional)')
    args = parser.parse_args()
    
    main(args.pickle_file, args.plot_output)
