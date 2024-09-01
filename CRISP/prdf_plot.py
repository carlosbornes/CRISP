import pickle
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation
from IPython.display import display, HTML

def plot_rdf_from_pickle(pickle_file):
    with open(pickle_file, 'rb') as f:
        data = pickle.load(f)
    
    x_data_all = data['x_data_all']
    y_data_all = data['y_data_all']

    # Calculate average RDF across all frames
    y_data_avg = np.mean(y_data_all, axis=0)

    # Find the index of the maximum y value in the average RDF
    max_y_index = np.argmax(y_data_avg)
    max_y_x = x_data_all[max_y_index]
    max_y = y_data_avg[max_y_index]

    # Plot the average RDF
    plt.plot(x_data_all, y_data_avg, label='Average RDF')
    plt.xlabel('Distance (Å)')
    plt.ylabel('Distribution')
    plt.title('Average Radial Distribution Function')

    # Add a vertical dashed line for the maximum y value
    plt.axvline(x=max_y_x, color='red', linestyle='--', label=f'Max Y ({max_y_x:.2f} Å)')

    # Add the legend
    plt.legend()
    plt.grid(True)
    plt.show()

def animate_rdf_from_pickle(pickle_file, save_file=None):
    with open(pickle_file, 'rb') as f:
        data = pickle.load(f)

    x_data_all = data['x_data_all']
    y_data_all = data['y_data_all']

    # Create a figure and axis
    fig, ax = plt.subplots()

    # Set x and y axis labels
    plt.xlabel('Distance (Å)')
    plt.ylabel('Distribution')

    def update(frame):
        ax.clear()  # Clear the previous frame
        y = y_data_all[frame]
        ax.plot(x_data_all, y)
        
        # Find the index of the maximum y value
        max_y_index = np.argmax(y)
        max_y_x = x_data_all[max_y_index]
        max_y = y[max_y_index]

        # Add a line for the maximum y value
        ax.axvline(x=max_y_x, color='red', linestyle='--', label=f'Max Y ({max_y_x:.2f} Å)')

        # Add grid lines
        ax.grid(True)

        # Update the legend with the frame number and max Y value info
        ax.legend([f'Frame {frame}', f'Max Y ({max_y_x:.2f} Å)'], loc='upper right')

    # Create the animation
    ani = FuncAnimation(fig, update, frames=range(len(y_data_all)), interval=500)

    if save_file:
        ani.save(save_file, writer='ffmpeg')

    # Display the animation in the Jupyter notebook
    display(HTML(ani.to_jshtml()))


