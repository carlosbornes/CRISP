import numpy as np
import ase.io
from ase.units import fs
import matplotlib.pyplot as plt
from ase.units import fs as fs_conversion
import statsmodels.api as sm


class DiffusionCoefficient:

    def __init__(self, traj, timestep, atom_indices=None, molecule=False):
        self.traj = traj
        self.timestep = timestep        

        # Condition used if user wants to calculate diffusion coefficients for specific atoms or all atoms
        self.atom_indices = atom_indices
        if self.atom_indices == None:
            self.atom_indices = [i for i in range(len(traj[0]))]

        # Condition if we are working with the mobility of a molecule, need to manage arrays slightly differently
        self.is_molecule = molecule
        if self.is_molecule:
            self.types_of_atoms = ["molecule"]
            self.no_of_atoms = [1]
        else:
            self.types_of_atoms = sorted(set(traj[0].symbols[self.atom_indices]))
            self.no_of_atoms = [len(self.atom_indices)]   
            
            
        # Dummy initialisation for important results data object
        self._slopes = []        
            
    @property
    def no_of_types_of_atoms(self):
        """

        Dynamically returns the number of different atoms in the system

        """
        return len(self.types_of_atoms)

    @property
    def slopes(self):
        """

        Method to return slopes fitted to datapoints. If undefined, calculate slopes

        """
        if len(self._slopes) == 0:
            self.calculate()
        return self._slopes

    @slopes.setter
    def slopes(self, values):
        """

        Method to set slopes as fitted to datapoints

        """
        self._slopes = values
        
        
    def _initialize_rmsd_dict(self):
        for atom_index in self.atom_indices:
            self.rmsd_values_by_atom[atom_index] = []

    def _initialise_arrays(self, ignore_n_images, number_of_segments):
        """

        Private function to initialise data storage objects. This includes objects to store the total timesteps
        sampled, the average diffusivity for species in any given segment, and objects to store gradient and intercept from fitting.

        Parameters:
            ignore_n_images (Int):
                Number of images you want to ignore from the start of the trajectory, e.g. during equilibration
            number_of_segments (Int):
                Divides the given trajectory in to segments to allow statistical analysis

        """

        total_images = len(self.traj) - ignore_n_images
        self.no_of_segments = number_of_segments
        self.len_segments = total_images // self.no_of_segments

        # These are the data objects we need when plotting information. First the x-axis, timesteps
        self.timesteps = np.linspace(0, total_images * self.timestep, total_images + 1)
        # This holds all the data points for the diffusion coefficients, averaged over atoms
        self.xyz_segment_ensemble_average = np.zeros(
            (self.no_of_segments, self.no_of_types_of_atoms, 3, self.len_segments))
        
        # This holds all the information on linear fits, from which we get the diffusion coefficients
        self.slopes = np.zeros((self.no_of_types_of_atoms, self.no_of_segments, 3))
        self.intercepts = np.zeros((self.no_of_types_of_atoms, self.no_of_segments, 3))

        self.cont_xyz_segment_ensemble_average = 0

    def calculate(self, ignore_n_images=0, number_of_segments=1):
        """

        Calculate the diffusion coefficients, using the previously supplied data. The user can break the data into segments and
        take the average over these trajectories, therefore allowing statistical analysis and derivation of standard deviations.
        Option is also provided to ignore initial images if data is perhaps unequilibrated initially.

        Parameters:
            ignore_n_images (Int):
                Number of images you want to ignore from the start of the trajectory, e.g. during equilibration
            number_of_segments (Int):
                Divides the given trajectory in to segments to allow statistical analysis

        """

        # Setup all the arrays we need to store information
        self._initialise_arrays(ignore_n_images, number_of_segments)

        for segment_no in range(self.no_of_segments):
            start = segment_no * self.len_segments
            end = start + self.len_segments
            seg = self.traj[ignore_n_images + start:ignore_n_images + end]
            
            # print("CHEK seg:::", len(seg))

            # If we are considering a molecular system, work out the COM for the starting structure
            if self.is_molecule:
                com_orig = np.zeros(3)
                for atom_no in self.atom_indices:
                    com_orig += seg[0].positions[atom_no] / len(self.atom_indices)

            for image_no in range(0, len(seg)):
                # This object collects the xyz displacements for all atom species in the image
                xyz_disp = np.zeros((self.no_of_types_of_atoms, 3))

                # Calculating for each atom individually, grouping by species type (e.g. solid state)
                if not self.is_molecule:
                    # For each atom, work out displacement from start coordinate and collect information with like atoms
                    for atom_no in self.atom_indices:                        
                        sym_index = self.types_of_atoms.index(seg[image_no].symbols[atom_no])                        
                        xyz_disp[sym_index] += np.square(seg[image_no].positions[atom_no] - seg[0].positions[atom_no])
                                                
                else:  # Calculating for group of atoms (molecule) and work out squared displacement
                    com_disp = np.zeros(3)
                    for atom_no in self.atom_indices:
                        com_disp += seg[image_no].positions[atom_no] / len(self.atom_indices)
                    xyz_disp[0] += np.square(com_disp - com_orig)
                
                
                # For each atom species or molecule, use xyz_disp to calculate the average data
                for sym_index in range(self.no_of_types_of_atoms):                 
                    denominator = (2 * self.no_of_atoms[sym_index])
                    
                    for xyz in range(3):
                        self.xyz_segment_ensemble_average[segment_no][sym_index][xyz][image_no] = (
                                    xyz_disp[sym_index][xyz] / denominator)

            # We've collected all the data for this entire segment, so now to fit the data.
            for sym_index in range(self.no_of_types_of_atoms):
                self.slopes[sym_index][segment_no], self.intercepts[sym_index][segment_no] = self._fit_data(
                    self.timesteps[start:end],
                    self.xyz_segment_ensemble_average[segment_no][sym_index])
                
        # Calculate RMSD for each particle at each time step
        self.rmsd_values = np.zeros((self.no_of_types_of_atoms, len(self.traj) - ignore_n_images))

        for image_no in range(len(self.traj)):
            for sym_index in range(self.no_of_types_of_atoms):
                initial_position = self.traj[ignore_n_images].positions[self.atom_indices[sym_index]]
                displacement = self.traj[image_no].positions[self.atom_indices[sym_index]] - initial_position
                squared_displacement = np.sum(displacement ** 2)
                rmsd = np.sqrt(squared_displacement / 3)  # Assuming 3D space
                self.rmsd_values[sym_index][image_no] = rmsd
                

    def _fit_data(self, x, y):
        """
        Private function that returns slope and intercept for linear fit to mean square diffusion data


        Parameters:
            x (Array of floats):
                Linear list of timesteps in the calculation
            y (Array of floats):
                Mean square displacement as a function of time.

        """

        # Initialise objects
        slopes = np.zeros(3)
        intercepts = np.zeros(3)

        # Convert into suitable format for lstsq
        x_edited = np.vstack([np.array(x), np.ones(len(x))]).T
        
        # Calculate slopes for x, y and z-axes
        for xyz in range(3):
            slopes[xyz], intercepts[xyz] = np.linalg.lstsq(x_edited, y[xyz], rcond=-1)[0]

        return slopes, intercepts
    

    def plot(self, ax=None, show=False):
        """

        Auto-plot of Diffusion Coefficient data. Provides basic framework for visualising analysis.

         Parameters:
            ax (Matplotlib.axes.Axes)
                Axes object on to which plot can be created
            show (Boolean)
                Whether or not to show the created plot. Default: False

        """

        # Create an x-axis that is in a more intuitive format for the view
        graph_timesteps = self.timesteps / fs_conversion

        for segment_no in range(self.no_of_segments):
            start = segment_no * self.len_segments
            end = start + self.len_segments
            label = None

                    
            ### total XYZ
            x = []
            y = []
            z = []
            
            for sym_index in range(self.no_of_types_of_atoms):
                for xyz in range(3):
                    if segment_no == 0:
                    # Add scatter graph  for the mean square displacement data in this segment  
                        x.append(self.xyz_segment_ensemble_average[segment_no][sym_index][0])
                        y.append(self.xyz_segment_ensemble_average[segment_no][sym_index][1])
                        z.append(self.xyz_segment_ensemble_average[segment_no][sym_index][2])
                    
        
               
        xt = graph_timesteps[:-1]
        y = (x[1]+y[1]+z[1])/3

        plt.loglog(xt, y, label="Original")
        plt.xlabel("Time (fs)")
        plt.ylabel("MSD_Fitted (Ang2)")
        plt.title("Log-Log Plot")
        plt.show()
        

        A = xt.reshape(-1, 1)  # Design matrix with only the x terms (no intercept)
        m, residuals, rank, s = np.linalg.lstsq(A, y, rcond=None)  # Fit y = mx

        diffusion_coefficient = m[0] 

        y_fitted = A @ m
        residual_std_error = np.sqrt(np.sum((y - y_fitted)**2) / (len(xt) - 1))  

        # Plot fitted line
        plt.plot(xt, y, label="Original MSD")
        plt.plot(xt, y_fitted, label="Fitted MSD")
        plt.xlabel("Time (fs)")
        plt.ylabel("Mean Square Displacement (Ang^2)")
        plt.title("MSD vs Time (Linear Fit)")
        plt.legend()
        plt.show()

        # Calculate standard error for diffusion coefficient (propagate from residuals)
        std_error = residual_std_error

        print("Diffusion Coefficient (m^2/sec^-1):", diffusion_coefficient * 10**(-5))
        print(f"Standard Error: {std_error * 10**(-5):.2e}")
        

        return diffusion_coefficient * 10**(-5), std_error * 10**(-5)
       


def main(args):
    traj = read(args.traj_file, index=args.traj_indices)
    diffusion_instance = DiffusionCoefficient(traj, args.timestep*fs, atom_indices)
    diffusion_instance.calculate()
    diffusion_instance.plot()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Calculate and plot Diffusion Coefficient from MD trajectory.")
    parser.add_argument("traj_file", type=str, help="Path to ASE-compatible trajectory file.")
    parser.add_argument("--atom_indices_file", type=str, help="Path to numpy array file containing atom indices.")
    parser.add_argument("--timestep", type=float, default=1.0, help="Timestep used in the simulation (default: 1.0 ps).")
    
    args = parser.parse_args()
    main(args)
