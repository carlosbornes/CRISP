# CRISP/data_analysis/msd.py

import ase.io
import numpy as np
import matplotlib.pyplot as plt
import csv
import argparse
import sys
from ase.units import fs
from ase.units import fs as fs_conversion
from scipy.optimize import curve_fit
import pandas as pd

def calculate_msd(traj, timestep, atom_indices=None, ignore_n_images=0):
    """
    Calculate Mean Square Displacement (MSD) vs time.
    
    Parameters:
        traj: Trajectory data
        timestep: Simulation timestep
        atom_indices: Indices of atoms to analyze (default: all atoms)
        ignore_n_images: Number of initial images to ignore
    
    Returns:
        msd_values: Mean square displacement values
        msd_times: Corresponding time values in femtoseconds
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
    
    Parameters:
        msd_times: Time values in femtoseconds
        msd_values: Mean square displacement values
        csv_file_path: Path to save the CSV file
    """
    with open(csv_file_path, 'w', newline='') as csvfile:
        csv_writer = csv.writer(csvfile)
        csv_writer.writerow(['Time (fs)', 'MSD'])
        for time, msd in zip(msd_times, msd_values):
            csv_writer.writerow([time, msd])
    
    print(f"MSD data has been saved to {csv_file_path}")

def calculate_diffusion_coefficient(msd_times, msd_values, start_index=None, end_index=None, 
                                   with_intercept=False, plot=False):
    """
    Calculate diffusion coefficient from MSD data.
    
    Parameters:
        msd_times: Time values in femtoseconds
        msd_values: Mean square displacement values
        start_index: Starting index for the fit
        end_index: Ending index for the fit
        with_intercept: Whether to fit with intercept
        plot: Whether to plot the fit
    
    Returns:
        D: Diffusion coefficient in cm²/s
        error: Statistical error in the diffusion coefficient
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
    if plot:
        plt.figure(figsize=(10, 6))
        plt.scatter(msd_times, msd_values, s=10, alpha=0.5, label='MSD data')
        plt.plot(x_fit, fit_func(x_fit), 'r-', linewidth=2, label=label)
        plt.xlabel('Time (fs)')
        plt.ylabel('MSD (Å²)')
        plt.title('Mean Square Displacement vs Time')
        plt.grid(True, alpha=0.3)
        plt.legend()
        plt.tight_layout()
        plt.show()
    
    print(f'Diffusion Coefficient: {D:.2e} cm²/s')
    print(f'Error: {error:.2e} cm²/s')
    
    return D, error

def calculate_save_msd(traj_file, timestep_value, indices_file=None, 
                      ignore_n_images=0, output_file="msd_results.csv", frame_skip=1):
    """
    Calculate MSD data and save to CSV file.
    
    Parameters:
        traj_file: Path to the ASE trajectory file
        timestep_value: Simulation timestep in ASE time units
        indices_file: Path to file containing atom indices
        ignore_n_images: Number of initial images to ignore
        output_file: Output CSV file path
        frame_skip: Number of frames to skip between samples
    
    Returns:
        msd_values: Mean square displacement values
        msd_times: Corresponding time values in femtoseconds
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

    # Adjust timestep for frame skipping
    timestep = timestep_value * frame_skip
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
                    with_intercept=False, plot=False):
    """
    Analyze MSD data from a CSV file.
    
    Parameters:
        csv_file: Path to the CSV file containing MSD data
        fit_start: Start index for fitting
        fit_end: End index for fitting
        with_intercept: Whether to fit with intercept
        plot: Whether to plot MSD vs time
    
    Returns:
        D: Diffusion coefficient
        error: Statistical error
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
            plot=plot
        )
        
        return D, error
        
    except Exception as e:
        print(f"Error analyzing MSD data: {e}")
        return None, None

def main():
    """Main function to run the MSD analysis from command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Calculate and analyze Mean Square Displacement (MSD) from trajectory data."
    )
    
    subparsers = parser.add_subparsers(dest="command", help="Command to run")
    
    # Calculate command
    calc_parser = subparsers.add_parser("calculate", help="Calculate and save MSD data")
    
    calc_parser.add_argument(
        "trajectory", 
        type=str, 
        help="Path to the ASE trajectory file."
    )
    
    calc_parser.add_argument(
        "--timestep", 
        type=float, 
        required=True,
        help="Simulation timestep in ASE time units (e.g., 0.5 * 100 * fs)"
    )
    
    calc_parser.add_argument(
        "--indices", 
        type=str, 
        help="Path to file containing atom indices (numpy array)",
        default=None
    )
    
    calc_parser.add_argument(
        "--ignore_images", 
        type=int, 
        default=0, 
        help="Number of initial images to ignore. Default is 0."
    )
    
    calc_parser.add_argument(
        "--output", 
        type=str, 
        default="msd_results.csv", 
        help="Output CSV file path. Default is 'msd_results.csv'."
    )
    
    calc_parser.add_argument(
        "--frame_skip", 
        type=int, 
        default=10, 
        help="Number of frames to skip between samples. Default is 10."
    )
    
    # Analyze command
    analyze_parser = subparsers.add_parser("analyze", help="Analyze saved MSD data")
    
    analyze_parser.add_argument(
        "--input", 
        type=str, 
        default="msd_results.csv", 
        help="Input CSV file path. Default is 'msd_results.csv'."
    )
    
    analyze_parser.add_argument(
        "--fit_start", 
        type=int,
        help="Start index for fitting"
    )
    
    analyze_parser.add_argument(
        "--fit_end", 
        type=int, 
        help="End index for fitting"
    )
    
    analyze_parser.add_argument(
        "--with_intercept", 
        action="store_true", 
        help="Fit with intercept"
    )
    
    analyze_parser.add_argument(
        "--plot", 
        action="store_true", 
        help="Plot MSD vs time"
    )
    
    args = parser.parse_args()
    
    if args.command == "calculate":
        # Convert timestep
        timestep = args.timestep * fs
        
        calculate_save_msd(
            traj_file=args.trajectory,
            timestep_value=timestep,
            indices_file=args.indices,
            ignore_n_images=args.ignore_images,
            output_file=args.output,
            frame_skip=args.frame_skip
        )
        
    elif args.command == "analyze":
        analyze_from_csv(
            csv_file=args.input,
            fit_start=args.fit_start,
            fit_end=args.fit_end,
            with_intercept=args.with_intercept,
            plot=args.plot
        )
        
    else:
        parser.print_help()


if __name__ == "__main__":
    main()