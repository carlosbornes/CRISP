import os
import numpy as np
from CRISP.data_analysis.h_bond import hydrogen_bonds, indices

# List of functionals in the specified order
functionals = ['M06-L', 'M11-L', 'MN12-L', 'MN15-L', 'revM06-L']

# Dictionary to map functional names to their directory names
func_dir_map = {
    'M06-L': 'M06L',
    'M11-L': 'M11L',
    'MN12-L': 'MN12L',
    'MN15-L': 'MN15L', 
    'revM06-L': 'REVM06L'
}

# Set parameters for hydrogen bond analysis
frame_skip = 10  # Skip every 10 frames
angle_cutoff = 120  # Minimum angle for hydrogen bond in degrees
h_bond_cutoff = 2.5  # Maximum H-bond distance in Angstroms
bond_cutoff = 1.6  # Maximum covalent bond distance in Angstroms
# time_step = 0.5  # Time step in fs (adjust as needed)

# Process each functional
for func in functionals:
    func_dir = func_dir_map[func]
    
    print(f"\nProcessing functional: {func}")
    
    # Define paths
    o_indices_path = f'./Data_supplementary_analysis/MetaGGA/{func_dir}/indices_new/O_indices.npy'
    h_indices_path = f'./Data_supplementary_analysis/MetaGGA/{func_dir}/indices_new/H_indices.npy'
    traj_path = f'./supplementary_data/MetaGGA/{func_dir}/TRAJEC.traj'
    output_dir = f'./Data_supplementary_analysis/MetaGGA/{func_dir}/hydrogen_bonds'
    
    os.makedirs(output_dir, exist_ok=True)
    
    if not os.path.exists(o_indices_path):
        print(f"Warning: Oxygen indices file not found for {func}: {o_indices_path}")
        continue
    
    if not os.path.exists(h_indices_path):
        print(f"Warning: Hydrogen indices file not found for {func}: {h_indices_path}")
        continue
    
    if not os.path.exists(traj_path):
        print(f"Warning: Trajectory file not found for {func}: {traj_path}")
        continue
    
    # Load oxygen and hydrogen indices
    try:
        o_indices = np.load(o_indices_path, allow_pickle=True)
        h_indices = np.load(h_indices_path, allow_pickle=True)
        print(f"Loaded {len(o_indices)} oxygen atoms and {len(h_indices)} hydrogen atoms")
    except Exception as e:
        print(f"Error loading indices for {func}: {e}")
        continue
    
    # Set up donor-acceptor indices path
    indices_path = f'{output_dir}/donor_acceptor_indices.npy'
    
    try:
        # Run hydrogen bond analysis
        print(f"Starting hydrogen bond analysis for {func}...")
        h_bonds_per_frame = hydrogen_bonds(
            traj_path=traj_path,
            frame_skip=frame_skip,
            acceptor_atoms=["O"],  # Use oxygen atoms as acceptors
            angle_cutoff=angle_cutoff,
            h_bond_cutoff=h_bond_cutoff,
            bond_cutoff=bond_cutoff,
            mic=True,
            single_h_bond=False,
            output_dir=output_dir,
            plot_count=True,
            plot_heatmap=True,
            plot_graph_frame=True,
            plot_graph_average=True,
            indices_path=indices_path,
            graph_frame_index=0
        )
        
        # Save hydrogen bonds per frame
        np.save(f'{output_dir}/h_bonds_per_frame.npy', np.array(h_bonds_per_frame))
        
        print(f"Completed hydrogen bond analysis for {func}")
        print(f"Average number of hydrogen bonds: {np.mean(h_bonds_per_frame):.2f}")
        print(f"Results saved to {output_dir}")
        