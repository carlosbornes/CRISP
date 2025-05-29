from CRISP.data_analysis.prdf import analyze_rdf
import numpy as np
import os

# Minnesota functionals from Villard et al. study
functionals = ['M06-L', 'M11-L', 'MN12-L', 'MN15-L', 'revM06-L']

# Directory mapping for file paths
func_dirs = {
    'M06-L': 'M06L',
    'M11-L': 'M11L', 
    'MN12-L': 'MN12L',
    'MN15-L': 'MN15L',
    'revM06-L': 'REVM06L'
}

peak_data = {}

print("Analyzing O-H radial distribution functions...")

for functional in functionals:
    dir_name = func_dirs[functional]
    print(f"\nWorking on {functional}...")
    
    # Set up file paths
    oxygen_indices = f'./Data_supplementary_analysis/MetaGGA/{dir_name}/indices_new/O_indices.npy'
    hydrogen_indices = f'./Data_supplementary_analysis/MetaGGA/{dir_name}/indices_new/H_indices.npy'
    trajectory = f'./supplementary_data/MetaGGA/{dir_name}/TRAJEC.traj'
    results_dir = f'./Data_supplementary_analysis/MetaGGA/{dir_name}/goh'
    
    os.makedirs(results_dir, exist_ok=True)
    
    try:
        # Load atom indices
        o_atoms = np.load(oxygen_indices)
        h_atoms = np.load(hydrogen_indices)
        
        # Setup for O-H partial RDF calculation
        atom_pairs = (o_atoms.tolist(), h_atoms.tolist())
        
        # Calculate RDF
        rdf_results = analyze_rdf(
            use_prdf=True,
            rmax=6.0,
            traj_path=trajectory,
            nbins=50,
            frame_skip=1,
            output_filename="prdf_o-atoms_h-atoms",
            atomic_indices=atom_pairs,
            output_dir=results_dir,
            create_plots=True
        )
        
        # Extract first coordination shell peak
        distances = rdf_results['x_data']
        avg_rdf = np.mean(rdf_results['y_data_all'], axis=0)
        first_peak_idx = np.argmax(avg_rdf)
        first_peak_pos = distances[first_peak_idx]
        
        peak_data[functional] = first_peak_pos
        print(f"  First O-H peak: {first_peak_pos:.2f} Å")
        
    except Exception as error:
        print(f"  Failed to process {functional}: {error}")

# Results summary
print("\n" + "="*40)
print("O-H COORDINATION ANALYSIS SUMMARY")
print("="*40)
print("Functional\t\tFirst Peak (Å)")
print("-"*40)

for func in functionals:
    if func in peak_data:
        tab = "\t\t" if len(func) <= 6 else "\t"
        print(f"{func}{tab}{peak_data[func]:.2f}")

print(f"\nAnalyzed {len(peak_data)} functionals successfully")