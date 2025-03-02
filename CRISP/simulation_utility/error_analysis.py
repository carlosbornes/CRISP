# CRISP/simulation_utility/error_analysis.py

import numpy as np
from statsmodels.tsa.stattools import acf
import warnings
import argparse

def optimal_lag(acf_values, threshold=0.05):
    """Determine optimal lag time based on autocorrelation threshold.
    
    Parameters
    ----------
    acf_values : numpy.ndarray
        Array of autocorrelation function values
    threshold : float, optional
        Threshold below which autocorrelation is considered negligible (default: 0.05)
        
    Returns
    -------
    int
        Optimal lag value where autocorrelation falls below threshold
    """
    for lag, value in enumerate(acf_values):
        if abs(value) < threshold:
            return lag
    
    acf_not_converged = (
    "Autocorrelation function is not converged. "
    f"Consider increasing the 'max_lag' parameter (current: {len(acf_values) - 1}) "
    "or extending the simulation length."
    )
    warnings.warn(acf_not_converged)
    
    return len(acf_values) - 1

def vector_acf(data, max_lag):
    """Calculate autocorrelation function for vector data.
    
    Parameters
    ----------
    data : numpy.ndarray
        2D array of vector data with shape (n_frames, n_dimensions)
    max_lag : int
        Maximum lag time for autocorrelation calculation
        
    Returns
    -------
    numpy.ndarray
        Array of autocorrelation values from lag 0 to max_lag
    """
    n_frames = data.shape[0]
    m = np.mean(data, axis=0)
    data_centered = data - m
    norm0 = np.mean(np.sum(data_centered**2, axis=1))
    acf_vals = np.zeros(max_lag + 1)
    for tau in range(max_lag + 1):
        dots = np.sum(data_centered[:n_frames - tau] * data_centered[tau:], axis=1)
        acf_vals[tau] = np.mean(dots) / norm0
    return acf_vals

def autocorrelation_analysis(data, max_lag=None, threshold=0.05):
    """Perform autocorrelation analysis on time series data.
    
    Parameters
    ----------
    data : numpy.ndarray
        1D array for scalar data or 2D array for vector data
    max_lag : int, optional
        Maximum lag time for autocorrelation calculation (default: min(1000, N/10))
    threshold : float, optional
        Threshold for determining optimal lag time (default: 0.05)
        
    Returns
    -------
    dict
        Dictionary containing:
        - mean: Mean value(s) of the data
        - acf_err: Statistical error estimated from autocorrelation
        - std: Standard deviation of the data
        - tau_int: Integrated autocorrelation time
        - optimal_lag: Optimal lag time determined by threshold
    """
    if data.ndim == 1:
        N = len(data)
        mean_value = np.mean(data)
        std_value = np.std(data, ddof=1)
        if max_lag is None:
            max_lag = min(1000, N // 10)
        acf_values = acf(data - mean_value, nlags=max_lag, fft=True)
        opt_lag = optimal_lag(acf_values, threshold)
    else:
        N = data.shape[0]
        mean_value = np.mean(data, axis=0)
        std_value = np.std(data, axis=0, ddof=1)
        if max_lag is None:
            max_lag = min(1000, N // 10)
        acf_values = vector_acf(data, max_lag)
        opt_lag = optimal_lag(acf_values, threshold)
    
    tau_int = 0.5 + np.sum(acf_values[1:opt_lag + 1])
    autocorr_error = std_value * np.sqrt(2 * tau_int / N)
    
    return {
        "mean": mean_value, 
        "acf_err": autocorr_error, 
        "std": std_value, 
        "tau_int": tau_int, 
        "optimal_lag": opt_lag
    }

def block_analysis(data, convergence_tol=0.001):
    """Perform block averaging analysis to estimate statistical errors.
    
    Parameters
    ----------
    data : numpy.ndarray
        1D array of time series data
    convergence_tol : float, optional
        Tolerance for determining convergence of block error (default: 0.001)
        
    Returns
    -------
    dict
        Dictionary containing:
        - mean: Mean value of the data
        - block_err: Statistical error estimated from block averaging
        - std: Standard deviation of the data
        - converged_blocks: Number of blocks at convergence
    """
    N = len(data)
    mean_value = np.mean(data)
    std_value = np.std(data, ddof=1)
    block_sizes = np.arange(1, N // 2)
    standard_errors = []
    
    for M in block_sizes:
        block_length = N // M
        truncated_data = data[:block_length * M]
        blocks = truncated_data.reshape(M, block_length)
        block_means = np.mean(blocks, axis=1)
        if len(block_means) > 1:
            std_error = np.std(block_means, ddof=1) / np.sqrt(M)
        else:
            continue
        
        standard_errors.append(std_error)
        if len(standard_errors) > 5:
            recent_errors = standard_errors[-5:]
            if np.max(recent_errors) - np.min(recent_errors) < convergence_tol:
                converged_blocks = M
                final_error = std_error
                break
    else:
        converged_blocks = block_sizes[-1]
        final_error = standard_errors[-1]
        warnings.warn("Block averaging did not fully converge. Consider increasing data length or lowering tolerance.")
    
    return {
        "mean": mean_value,
        "block_err": final_error,
        "std": std_value,
        "converged_blocks": converged_blocks
    }

def main():
    """Run autocorrelation and block analysis on input data file."""
    parser = argparse.ArgumentParser()
    parser.add_argument("input_file", type=str)
    parser.add_argument("--max_lag", type=int, default=500)
    parser.add_argument("--threshold", type=float, default=0.05)
    parser.add_argument("--convergence_tol", type=float, default=0.001)
    args = parser.parse_args()
    
    data_energy = np.loadtxt(args.input_file, skiprows=1, usecols=2)
    acf_error = autocorrelation_analysis(data_energy, max_lag=args.max_lag, threshold=args.threshold)
    block_error = block_analysis(data_energy, convergence_tol=args.convergence_tol)
    
    # Print results to console
    print("\nAutocorrelation Analysis Results:")
    print(f"Mean: {acf_error['mean']:.6f}")
    print(f"Statistical Error (ACF): {acf_error['acf_err']:.6e}")
    print(f"Integrated Correlation Time: {acf_error['tau_int']:.2f}")
    print(f"Optimal Lag: {acf_error['optimal_lag']}")
    
    print("\nBlock Analysis Results:")
    print(f"Mean: {block_error['mean']:.6f}")
    print(f"Statistical Error (Block): {block_error['block_err']:.6e}")
    print(f"Converged at {block_error['converged_blocks']} blocks")

if __name__ == "__main__":
    main()
