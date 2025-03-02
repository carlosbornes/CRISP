# CRISP/simulation_utility/atomic_traj_linemap.py

import os
import numpy as np
import argparse
from typing import List, Optional, Dict
from ase.io import read
import plotly.graph_objects as go

# Dictionary of covalent radii for common elements in Ångström
# Used to scale atom sizes appropriately in the visualization
COVALENT_RADII = {
    'H': 0.31, 'He': 0.28, 'Li': 1.28, 'Be': 0.96, 'B': 0.84, 'C': 0.76, 
    'N': 0.71, 'O': 0.66, 'F': 0.57, 'Ne': 0.58, 'Na': 1.66, 'Mg': 1.41, 
    'Al': 1.21, 'Si': 1.11, 'P': 1.07, 'S': 1.05, 'Cl': 1.02, 'Ar': 1.06, 
    'K': 2.03, 'Ca': 1.76, 'Fe': 1.24, 'Ni': 1.21, 'Cu': 1.38, 'Zn': 1.31
}

# Default color palette for atom types
ELEMENT_COLORS = {
    'H': 'white', 'C': 'black', 'N': 'blue', 'O': 'red', 'F': 'green',
    'Na': 'purple', 'Mg': 'pink', 'Al': 'gray', 'Si': 'yellow', 'P': 'orange',
    'S': 'yellow', 'Cl': 'green', 'K': 'purple', 'Ca': 'gray', 'Fe': 'orange',
    'Cu': 'orange', 'Zn': 'gray'
}

def plot_atomic_trajectory(
    traj_path: str,
    selected_indices: List[int],
    output_path: str,
    frame_skip: int = 100,
    plot_title: str = None,
    show_plot: bool = False,
    atom_size_scale: float = 1.0
):
    """
    Create a 3D visualization of atom trajectories with all atom types displayed.
    
    Parameters
    ----------

    traj_path : str
        Path to the ASE trajectory file
    selected_indices : List[int]
        Atom indices to plot trajectories for
    output_path : str
        Path to save the generated HTML visualization
    frame_skip : int, optional
        Use every nth frame from the trajectory (default: 100)
    plot_title : str, optional
        Custom title for the plot (default: auto-generated)
    show_plot : bool, optional
        Whether to display the plot (default: False)
    atom_size_scale : float, optional
        Scale factor for atom sizes (default: 1.0)
    """
    # Load trajectory
    print(f"Loading trajectory from {traj_path} (using every {frame_skip}th frame)...")
    traj = read(traj_path, index=f'::{frame_skip}')
    
    # Convert to list if not already (happens with single frame)
    if not isinstance(traj, list):
        traj = [traj]
    
    print(f"Loaded {len(traj)} frames from trajectory")
    
    if not selected_indices:
        print("No atom indices selected for trajectory plotting")
    else:
        print(f"Selected {len(selected_indices)} atoms for trajectory plotting: {selected_indices}")
    
    # Get box dimensions
    box = traj[0].cell.lengths()
    print(f"Simulation box dimensions: {box} Å")
    
    # Find all unique atom types in the first frame
    atom_types = {}
    max_index = max([atom.index for atom in traj[0]])
    print(f"Analyzing atom types in first frame (total atoms: {len(traj[0])}, max index: {max_index})...")
    
    for atom in traj[0]:
        symbol = atom.symbol
        if symbol not in atom_types:
            atom_types[symbol] = []
        atom_types[symbol].append(atom.index)
    
    print(f"Found {len(atom_types)} atom types: {', '.join(atom_types.keys())}")
    
    # Create figure
    fig = go.Figure()
    
    # Determine colors for selected indices
    # Use same color if more than 5 indices
    use_same_color = len(selected_indices) > 5
    colors = ['blue'] * len(selected_indices) if use_same_color else [
        'blue', 'green', 'red', 'orange', 'purple'
    ][:len(selected_indices)]
    
    # Add atoms of each type from the first frame
    for symbol, indices in atom_types.items():
        positions = np.array([traj[0].positions[i] for i in indices])
        
        # Skip if no atoms of this type
        if len(positions) == 0:
            continue
        
        # Determine marker size based on covalent radius
        size = COVALENT_RADII.get(symbol, 1.0) * 3.0 * atom_size_scale
        color = ELEMENT_COLORS.get(symbol, 'gray')
        
        fig.add_trace(go.Scatter3d(
            x=positions[:, 0],
            y=positions[:, 1],
            z=positions[:, 2],
            mode='markers',
            name=f'{symbol} Atoms',
            marker=dict(
                size=size,
                color=color,
                symbol='circle',
                opacity=0.7,
                line=dict(color='black', width=0.5)
            )
        ))
    
    # Collect trajectories for selected atoms
    selected_positions = {idx: [] for idx in selected_indices}
    for atoms in traj:
        for idx in selected_indices:
            if idx < len(atoms):
                selected_positions[idx].append(atoms.positions[idx])
            else:
                print(f"Warning: Index {idx} is out of range")
    
    # Add trajectories for selected atoms
    annotations = []
    for i, idx in enumerate(selected_indices):
        if not selected_positions[idx]:  # Skip if no positions collected
            continue
            
        pos = np.array(selected_positions[idx])
        color = colors[i % len(colors)]
        
        # Add trajectory line and markers
        fig.add_trace(go.Scatter3d(
            x=pos[:, 0],
            y=pos[:, 1],
            z=pos[:, 2],
            mode='lines+markers',
            name=f'Atom {idx}',
            line=dict(width=3, color=color),
            marker=dict(size=4, color=color),
        ))
        
        # Add annotations for first and last frames
        first_frame_pos = pos[0]
        last_frame_pos = pos[-1]
        
        annotations.extend([
            dict(
                x=first_frame_pos[0],
                y=first_frame_pos[1],
                z=first_frame_pos[2],
                text=f'Start {idx}',
                showarrow=True,
                arrowhead=2,
                ax=20,
                ay=-20,
                arrowcolor=color,
                font=dict(color=color, size=10)
            ),
            dict(
                x=last_frame_pos[0],
                y=last_frame_pos[1],
                z=last_frame_pos[2],
                text=f'End {idx}',
                showarrow=True,
                arrowhead=2,
                ax=-20,
                ay=-20,
                arrowcolor=color,
                font=dict(color=color, size=10)
            )
        ])
    
    # Set plot title
    if not plot_title:
        atom_types_str = ', '.join(atom_types.keys())
        plot_title = f'Atomic Trajectories in {atom_types_str} System'
    
    # Update layout
    fig.update_layout(
        title=plot_title,
        scene=dict(
            xaxis_title='X (Å)',
            yaxis_title='Y (Å)',
            zaxis_title='Z (Å)',
            xaxis=dict(range=[0, box[0]]),
            yaxis=dict(range=[0, box[1]]),
            zaxis=dict(range=[0, box[2]]),
            aspectmode='cube'
        ),
        margin=dict(l=0, r=0, b=0, t=40),
        scene_annotations=annotations
    )
    
    # Create output directory if needed
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    
    # Save figure
    fig.write_html(output_path)
    print(f"Plot has been saved to {output_path}")
    
    # Show figure if requested
    if show_plot:
        fig.show()
    
    return fig

def parse_indices(indices_str):
    """Parse comma-separated indices."""
    if not indices_str:
        return []
    
    try:
        return [int(idx.strip()) for idx in indices_str.split(',')]
    except ValueError:
        raise argparse.ArgumentTypeError("Indices must be comma-separated integers")

def main():
    """Parse command-line arguments and create the trajectory visualization."""
    parser = argparse.ArgumentParser()
    
    # Required arguments
    parser.add_argument("traj_path", type=str)
    parser.add_argument("output_path", type=str)
    
    # Optional arguments
    parser.add_argument("--selected-indices", type=parse_indices, default=[])
    parser.add_argument("--frame-skip", type=int, default=100)
    parser.add_argument("--plot-title", type=str, default=None)
    parser.add_argument("--show-plot", action="store_true")
    parser.add_argument("--atom-size-scale", type=float, default=1.0)
    
    args = parser.parse_args()
    
    # Ensure output path ends with .html
    if not args.output_path.endswith('.html'):
        args.output_path += '.html'
    
    # Create the plot
    plot_atomic_trajectory(
        args.traj_path,
        args.selected_indices,
        args.output_path,
        args.frame_skip,
        args.plot_title,
        args.show_plot,
        args.atom_size_scale
    )

if __name__ == "__main__":
    main()
