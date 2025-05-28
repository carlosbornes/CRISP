from error_analysis import error_analysis
import numpy as np

# MD trajectory energy
data = np.load("energy.npy")
res = error_analysis(data)
print(res["acf_results"])
print(res["block_results"])

# Predicted chemical shift for each frame
data = np.loadtxt("Al_nmr_results.txt",usecols=3,skiprows=1)
res = error_analysis(data)
print(res["acf_results"])
print(res["block_results"])

# Error of the Al-O-Si angles
data = np.loadtxt("Al_nmr_results.txt",usecols=9,skiprows=1)
res = error_analysis(data)
print(res["acf_results"])
print(res["block_results"])
