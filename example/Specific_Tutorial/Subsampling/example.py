from subsampling import subsample

# Filename of the trajectory which to subsample
# Can be either individual one or *.xyz or *.traj to combine every file
filename = "rMD17.xyz"

# Number of subsampled structures (Default: 50)
# Has to be lower or equal to number of structures
n_samples = 10

# Atomic indices to focus the subsampling (Default: "all")
# Options
### all atoms: "all"
### element types: "O" or ["O","H"] or ("O","H")
### specific indices: 1 or [1,2,3]
### NumPy file with an array: "indices.npy"
index_type = "all"

# File format of the trajectories (Default: None)
file_format = "extxyz"

# How many frames to skip (Default: 1)
skip = 1

# Plot the convergence of the subsampling (Default: False)
plot_subsample = True

###

# Example use (outputs subsampled trajectory)
subsample(filename=filename,
          n_samples=n_samples,
          index_type=index_type,
          file_format=file_format,
          skip=skip,
          plot_subsample=plot_subsample)

###


