"""
CRISP/data_analysis/msd.py

This module performs mean square displacement (MSD) analysis on molecular dynamics 
trajectory data for diffusion coefficient calculations.
"""

import ase.io
import numpy as np
import matplotlib.pyplot as plt
import csv
import sys
from ase.units import fs
from ase.units import fs as fs_conversion
from scipy.optimize import curve_fit
import pandas as pd
import os

def calculate_msd(traj, timestep, atom_indices=None, ignore_n_images=0):
    """
    Calculate Mean Square Displacement (MSD) vs time.
    
    Parameters
    ----------
    traj : list of ase.Atoms
        Trajectory data
    timestep : float
        Simulation timestep
    atom_indices : numpy.ndarray, optional
        Indices of atoms to analyze (default: all atoms)
    ignore_n_images : int, optional
        Number of initial images to ignore (default: 0)
    
    Returns
    -------
    tuple
        (msd_values, msd_times) where msd_values are Mean square displacement values
        and msd_times are corresponding time values in femtoseconds
    """
    # Initialize parameters
    if atom_indices is None:
        atom_indices = list(range(len(traj[0])))
    
    # Calculate time values
    total_images = len(traj) - ignore_n_images
    timesteps = np.linspace(0, total_images * timestep, total_images + 1)
    msd_times = timesteps[:-1] / fs_conversion  # Convert to femtoseconds
    
    # Calculate MSD
    msd_values = np.zeros(total_images)
    for i in range(ignore_n_images, len(traj)):
        idx = i - ignore_n_images
        displacements = np.zeros(3)
        
        for atom_idx in atom_indices:
            displacement = traj[i].positions[atom_idx] - traj[ignore_n_images].positions[atom_idx]
            displacements += np.square(displacement)
            
        # Average MSD over all atoms
        msd_values[idx] = np.sum(displacements) / (3 * len(atom_indices))
    
    return msd_values, msd_times

def save_msd_data(msd_times, msd_values, csv_file_path):
    """
    Save MSD data to a CSV file.
    
    Parameters
    ----------
    msd_times : numpy.ndarray
        Time values in femtoseconds
    msd_values : numpy.ndarray
        Mean square displacement values
    csv_file_path : str
        Path to save the CSV file
        
    Returns
    -------
    None
    """
    with open(csv_file_path, 'w', newline='') as csvfile:
        csv_writer = csv.writer(csvfile)
        csv_writer.writerow(['Time (fs)', 'MSD'])
        for time, msd in zip(msd_times, msd_values):
            csv_writer.writerow([time, msd])
    
    print(f"MSD data has been saved to {csv_file_path}")

def calculate_diffusion_coefficient(msd_times, msd_values, start_index=None, end_index=None, 
                                   with_intercept=False, plot_msd=False):
    """
    Calculate diffusion coefficient from MSD data.
    
    Parameters
    ----------
    msd_times : numpy.ndarray
        Time values in femtoseconds
    msd_values : numpy.ndarray
        Mean square displacement values
    start_index : int, optional
        Starting index for the fit (default: 1/3 of data length)
    end_index : int, optional
        Ending index for the fit (default: None)
    with_intercept : bool, optional
        Whether to fit with intercept (default: False)
    plot_msd : bool, optional
        Whether to plot the fit (default: False)
    
    Returns
    -------
    tuple
        (D, error) where D is the diffusion coefficient in cm²/s and error is the statistical error
    """
    # Set default indices if not provided
    if start_index is None:
        start_index = len(msd_times) // 3  # Start at 1/3 of the data by default
    
    if end_index is None:
        end_index = len(msd_times)
    
    # Ensure the indices are valid
    if start_index < 0 or end_index > len(msd_times):
        raise ValueError("Indices are out of bounds.")
    if start_index >= end_index:
        raise ValueError("Start index must be less than end index.")
    
    # Get the data to fit
    x_fit = msd_times[start_index:end_index]
    y_fit = msd_values[start_index:end_index]
    
    # Define the linear models
    def linear_no_intercept(x, m):
        return m * x
    
    def linear_with_intercept(x, m, c):
        return m * x + c
    
    # Perform the fit
    if with_intercept:
        params, covariance = curve_fit(linear_with_intercept, x_fit, y_fit)
        slope, intercept = params
        fit_func = lambda x: linear_with_intercept(x, slope, intercept)
        label = f'D = {slope * 10**-5:.2e} cm²/s'
    else:
        params, covariance = curve_fit(linear_no_intercept, x_fit, y_fit)
        slope = params[0]
        intercept = 0
        fit_func = lambda x: linear_no_intercept(x, slope)
        label = f'D = {slope * 10**-5:.2e} cm²/s'
    
    # Calculate error
    std_err = np.sqrt(np.diag(covariance))[0]
    
    # Calculate diffusion coefficient (D = slope/6 * 10^-5)
    D = slope * 10**-5
    error = std_err * 10**-5
    
    # Plot if requested
    if plot_msd:
        plt.figure(figsize=(10, 6))
        plt.scatter(msd_times, msd_values, s=10, alpha=0.5, label='MSD data')
        plt.plot(x_fit, fit_func(x_fit), 'r-', linewidth=2, label=label)
        plt.xlabel('Time (fs)')
        plt.ylabel('MSD (Å²)')
        plt.title('Mean Square Displacement vs Time')
        plt.grid(True, alpha=0.3)
        plt.legend()
        plt.tight_layout()
        plt.savefig('msd_analysis.png', dpi=300, bbox_inches='tight')
        plt.show()
    
    print(f'Diffusion Coefficient: {D:.2e} cm²/s')
    print(f'Error: {error:.2e} cm²/s')
    
    return D, error

def calculate_save_msd(traj_file, timestep_value, indices_file=None, 
                      ignore_n_images=0, output_file="msd_results.csv", frame_skip=1):
    """
    Calculate MSD data and save to CSV file.
    
    Parameters
    ----------
    traj_file : str
        Path to the ASE trajectory file
    timestep_value : float
        Simulation timestep in ASE time units
    indices_file : str, optional
        Path to file containing atom indices (default: None)
    ignore_n_images : int, optional
        Number of initial images to ignore (default: 0)
    output_file : str, optional
        Output CSV file path (default: "msd_results.csv")
    frame_skip : int, optional
        Number of frames to skip between samples (default: 1)
    
    Returns
    -------
    tuple
        (msd_values, msd_times) - MSD values and corresponding time values
    """
    # Load trajectory file
    try:
        full_traj = ase.io.read(traj_file, ':')
        print(f"Loaded full trajectory with {len(full_traj)} frames")
        
        # Create new trajectory with frame skipping
        traj = []
        for i in range(0, len(full_traj), frame_skip):
            traj.append(full_traj[i])
        
        print(f"Using {len(traj)} frames after skipping every {frame_skip} frames")
    except Exception as e:
        print(f"Error loading trajectory file: {e}")
        return None, None

    # Load atom indices if specified
    atom_indices = None
    if indices_file:
        try:
            atom_indices = np.load(indices_file)
            print(f"Loaded {len(atom_indices)} atom indices")
        except Exception as e:
            print(f"Error loading atom indices: {e}")
            return None, None

    timestep = timestep_value 
    print(f"Using adjusted timestep: {timestep/fs} * fs (original: {timestep_value/fs} * fs)")

    # Calculate MSD
    print("Calculating MSD...")
    msd_values, msd_times = calculate_msd(
        traj=traj, 
        timestep=timestep,
        atom_indices=atom_indices,
        ignore_n_images=ignore_n_images
    )

    # Save MSD data
    save_msd_data(
        msd_times=msd_times, 
        msd_values=msd_values, 
        csv_file_path=output_file
    )
    
    return msd_values, msd_times

def analyze_from_csv(csv_file="msd_results.csv", fit_start=None, fit_end=None, 
                    with_intercept=False, plot_msd=False):
    """
    Analyze MSD data from a CSV file.
    
    Parameters
    ----------
    csv_file : str, optional
        Path to the CSV file containing MSD data (default: "msd_results.csv")
    fit_start : int, optional
        Start index for fitting (default: None)
    fit_end : int, optional
        End index for fitting (default: None)
    with_intercept : bool, optional
        Whether to fit with intercept (default: False)
    plot_msd : bool, optional
        Whether to plot MSD vs time (default: False)
    
    Returns
    -------
    tuple
        (D, error) where D is the diffusion coefficient in cm²/s and error is the statistical error
    """
    try:
        # Load MSD data from CSV
        df = pd.read_csv(csv_file)
        print(f"Loaded MSD data from {csv_file}")
        
        # Extract time and MSD values
        msd_times = df['Time (fs)'].values
        msd_values = df['MSD'].values
        
        # Calculate diffusion coefficient
        D, error = calculate_diffusion_coefficient(
            msd_times=msd_times, 
            msd_values=msd_values, 
            start_index=fit_start, 
            end_index=fit_end, 
            with_intercept=with_intercept, 
            plot_msd=plot_msd
        )
        
        return D, error
        
    except Exception as e:
        print(f"Error analyzing MSD data: {e}")
        return None, None

def msd_analysis(traj_file, timestep_value, indices_file=None, ignore_n_images=0,
                output_dir="msd_results", frame_skip=10, fit_start=None, fit_end=None,
                with_intercept=False, plot_msd=True):
    """
    Perform complete MSD analysis workflow: calculate MSD, save data, and analyze diffusion.
    
    Parameters
    ----------
    traj_file : str
        Path to the ASE trajectory file
    timestep_value : float
        Simulation timestep in ASE time units
    indices_file : str, optional
        Path to file containing atom indices (default: None)
    ignore_n_images : int, optional
        Number of initial images to ignore (default: 0)
    output_dir : str, optional
        Directory to save output files (default: "msd_results")
    frame_skip : int, optional
        Number of frames to skip between samples (default: 10)
    fit_start : int, optional
        Start index for fitting diffusion coefficient (default: None)
    fit_end : int, optional
        End index for fitting diffusion coefficient (default: None)
    with_intercept : bool, optional
        Whether to fit with intercept (default: False)
    plot_msd : bool, optional
        Whether to plot results (default: True)
        
    Returns
    -------
    dict
        Dictionary containing MSD values, times, diffusion coefficient and error
    """
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Define output CSV path
    csv_path = os.path.join(output_dir, "msd_data.csv")
    
    # Convert timestep to ASE units if needed
    if timestep_value < 1:  # Assuming user provided in fs
        timestep = timestep_value * fs
    else:
        timestep = timestep_value
    
    # Calculate and save MSD data
    msd_values, msd_times = calculate_save_msd(
        traj_file=traj_file,
        timestep_value=timestep,
        indices_file=indices_file,
        ignore_n_images=ignore_n_images,
        output_file=csv_path,
        frame_skip=frame_skip
    )
    
    if msd_values is None or msd_times is None:
        print("Failed to calculate MSD.")
        return None
    
    # Calculate diffusion coefficient
    D, error = calculate_diffusion_coefficient(
        msd_times=msd_times,
        msd_values=msd_values,
        start_index=fit_start,
        end_index=fit_end,
        with_intercept=with_intercept,
        plot_msd=plot_msd
    )
    
    # Save summary to text file
    summary_path = os.path.join(output_dir, "msd_summary.txt")
    with open(summary_path, 'w') as f:
        f.write("Mean Square Displacement (MSD) Analysis Summary\n")
        f.write("=============================================\n\n")
        f.write(f"Trajectory file: {traj_file}\n")
        f.write(f"Timestep: {timestep/fs} fs\n")
        f.write(f"Frame skip: {frame_skip}\n")
        if indices_file:
            f.write(f"Atom indices file: {indices_file}\n")
        f.write(f"Number of frames analyzed: {len(msd_values)}\n\n")
        f.write(f"Diffusion Coefficient: {D:.6e} cm²/s\n")
        f.write(f"Statistical Error: {error:.6e} cm²/s\n")
        if fit_start is not None and fit_end is not None:
            f.write(f"Fit range: {fit_start} to {fit_end}\n")
        f.write(f"Fit with intercept: {with_intercept}\n")
    
    print(f"Analysis summary saved to {summary_path}")
    
    # Return results as dictionary
    return {
        "msd_values": msd_values,
        "msd_times": msd_times,
        "diffusion_coefficient": D,
        "error": error
    }