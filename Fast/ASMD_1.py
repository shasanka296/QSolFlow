import subprocess
import re
import os
import gromacs
import MDAnalysis as mda
import MDAnalysis.analysis.rdf as rdf
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import pandas as pd
from d3tales_fw.workflows.envwf import meta_dir, MDP_Location


class ASMD:

    """A class that uses the GROMACS wrapper along with subprocess to create and run simulation using
    GROMACS 2020 software package to find energies and uses pandas,
    NumPy and matplotlib to store and visualize data, as well as to generate  RDF plots.
    A density check function is built into the class though, it is disabled due to
    lack of reliable data about solution densities for arbitrary mixtures.
    The order of function class for the standard simulation is as follows:
    count_warnings,
    EnergyMin,
    NVT,
    NPT,
    calculate_density,
    check_density_accuracy (NOTE: CURRENTLY THE FUNCTION IS CALLED FOR THE SAKE OF CONSISTENCY, IT DOES NOT ACTUALY VERIFY THE DENSITY),
    correct,
    index_file,
    extract_residues_from_itp,
    RDF,
    coordination_number.
    """



    def __init__(self,key):

        self.current_dir = meta_dir
        self.input_dir = os.path.join(self.current_dir,f"InputGrofiles{key}")
        self.output_dir =  os.path.join(self.current_dir,f"Output{key}")
        self.mdp_dir = os.path.join(self.current_dir,f"InputGrofiles{key}","MDP")
        self.topology_file = os.path.join(self.input_dir,"topol.top")
        self.initial_coordinates =os.path.join(self.input_dir,"solvated.gro")
        self.xtc_file = os.path.join(self.output_dir,f"production.xtc")
        self.tpr_file = os.path.join(self.output_dir,f"production.tpr")
        self.eqtpr_file = os.path.join(self.output_dir,f"equilibration.tpr")
        self.trr_file = os.path.join(self.output_dir,f"production.trr")
        self.nmol_itp = os.path.join(self.input_dir,"nmol.itp")
        self.gro_file =os.path.join(self.output_dir,f"production.gro")
        self.edr_file = os.path.join(self.output_dir,f"equilibration.edr")


        # Copy MDP files to output directory
        subprocess.run(["cp", "-r", f"{self.mdp_dir}/.", self.output_dir])
        # Create output directory if it doesn't exist
        os.makedirs(self.output_dir, exist_ok=True)
        # Set initial values
        self.max_warn = 10
        self.threshold = 100

    # Run simulations

#1
    def count_warnings(self,filename):
        # Read the simulation output log file and count the number of warnings
        if os.path.isfile(self.output_dir + "/em.log") == True:
            with open(filename, 'r') as file:
                log_data = file.read()

            warning_pattern = r'WARNING'
            warning_count = len(re.findall(warning_pattern, log_data))
        else:
            warning_count = 100

        return warning_count


    def update_maxwarn(self,max_warn, threshold):
        # Check the warning count and update the max_warn value if needed
        warning_count = self.count_warnings(os.path.join(self.output_dir, "em.log"))
        if warning_count > threshold:
            max_warn *= 2  # Double the max_warn value
            print(f"Warning count ({warning_count}) exceeded threshold. Updating max_warn to: {max_warn}")
        else:
            print(f"Simulation completed successfully with {warning_count} warnings.")

        return max_warn


    def get_available_cpu_threads(self):
        # Get the number of available CPU threads
        cpu_info = os.popen('lscpu').read()
        threads_pattern = r'Thread\(s\) per core:\s+(\d+)'
        result = re.search(threads_pattern, cpu_info)
        if result:
            threads_per_core = int(result.group(1))
            available_threads = os.cpu_count() * threads_per_core
            return available_threads
        else:
            return None


    def get_available_gpus(self):
        # Get the number of available GPUs
        gpu_info = os.popen('nvidia-smi --list-gpus').read()
        gpu_count = len(re.findall(r'GPU\s\d:', gpu_info))
        return gpu_count


    def run_gromacs_simulation(self,command, max_warn, threshold):
        # Execute the GROMACS simulation command
        subprocess.run(command + ["-maxwarn", str(max_warn)])

        # Check the warning count and update max_warn if needed
        while True:
            max_warn = self.update_maxwarn(max_warn, threshold)
            if self.count_warnings(os.path.join(self.output_dir, "em.log")) <= threshold:
                break

        return max_warn
#2
    def EnergyMin(self):
        command = ["gmx_mpi", "grompp", "-f", os.path.join(self.mdp_dir, "em.mdp"), "-c", self.initial_coordinates, "-p", self.topology_file,
                   "-o", os.path.join(self.output_dir, "em.tpr")]
        self.run_gromacs_simulation(command, self.max_warn, self.threshold)

        available_threads = self.get_available_cpu_threads()
        available_gpus = self.get_available_gpus()

        command = ["gmx_mpi", "mdrun", "-deffnm", "em", "-v"]
        if available_threads is not None:
            command += ["-ntomp", str(available_threads)]
        subprocess.run(command, cwd=self.output_dir)
#3
    def NVT(self):
        command = ["gmx_mpi", "grompp", "-f", os.path.join(self.mdp_dir, "nvt.mdp"), "-c", os.path.join(self.output_dir, "em.gro"), "-p",
                   self.topology_file, "-o", os.path.join(self.output_dir, "nvt.tpr")]
        self.run_gromacs_simulation(command, self.max_warn, self.threshold)

        command = ["gmx_mpi", "mdrun", "-deffnm", "nvt", "-v"]
        if self.get_available_cpu_threads() is not None:
            command += ["-ntomp", str(self.get_available_cpu_threads())]
        subprocess.run(command, cwd=self.output_dir)
#4
    def NPT(self):
        command = ["gmx_mpi", "grompp", "-f", os.path.join(self.mdp_dir, "equilibration.mdp"), "-c",
                   os.path.join(self.output_dir, "nvt.gro"), "-t", os.path.join(self.output_dir, "nvt.cpt"), "-p", self.topology_file, "-o",
                   os.path.join(self.output_dir, "equilibration.tpr")]
        self.run_gromacs_simulation(command, self.max_warn, self.threshold)
        command = ["gmx_mpi", "mdrun", "-deffnm", "equilibration", "-v"]
        # if available_gpus > 0:
        # command = ["mpirun", "-np", str(available_gpus), "gmx_mpi", "mdrun", "-deffnm", "equilibration", "-v", "-nb", "gpu"]
        if self.get_available_cpu_threads() is not None:
            command += ["-ntomp", str(self.get_available_cpu_threads())]
        subprocess.run(command, cwd=self.output_dir)

    # Calculate density for validation


#5
    def calculate_density(self):
        # Ensure GROMACS commands are available.
        gromacs.config.setup()

        # Run gmx_mpi energy to extract density.
        density_xvg_file = f"{self.output_dir}/density.xvg"
        # command = f'gmx_mpi energy -s {os.path.abspath(eqtpr_file)} -f {os.path.abspath(edr_file)} -o {output_dir}/{density_xvg_file}'
        energy_mpi = gromacs.tools.Energy_mpi(s=os.path.abspath(self.eqtpr_file),
                                              f=os.path.abspath(self.edr_file),
                                              o=f'{density_xvg_file}')

        energy_mpi.run(input="Density")
        # subprocess.run(command)

        densities = []
        with open(f"{self.output_dir}/density.xvg", "r") as f:
            for line in f.readlines():
                if line.startswith("#") or line.startswith("@"):
                    continue
                densities.append(float(line.split()[1]))

        return densities

# rubn this as is
    def extract_density(self):
        # Read density values from the xvg file.
        densities = []
        with open(f"{self.output_dir}/density.xvg", "r") as f:
            for line in f.readlines():
                if line.startswith("#") or line.startswith("@"):
                    continue
                densities.append(float(line.split()[1]))

        return densities

    def utlizie(self):
        self.densities = self.extract_density()

        print("Density values:", self.densities)
        return self.densities
       # Run production if density is within accuracy
    def production_run(self,topology_file, output_dir):
        max_warn = 10
        threshold = 100
        command = ["gmx_mpi", "grompp", "-f", os.path.join(self.mdp_dir, "production.mdp"), "-c",
                   os.path.join(output_dir, "equilibration.gro"), "-t", os.path.join(output_dir, "equilibration.cpt"), "-p",
                   topology_file, "-o",
                   os.path.join(output_dir, "production.tpr")]
        self.run_gromacs_simulation(command, max_warn, threshold)
        command = ["gmx_mpi", "mdrun", "-deffnm", "production", "-v"]
        # if available_gpus > 0:
        # command = ["mpirun", "-np", str(available_gpus), "gmx_mpi", "mdrun", "-deffnm", "production", "-v", "-nb", "gpu"]
        if self.get_available_cpu_threads() is not None:
            command += ["-ntomp", str(self.get_available_cpu_threads())]
        subprocess.run(command, cwd= self.output_dir)
        print("Simulation completed.")


    def check_density_accuracy(self,x, densities):
        mean_density = sum(densities) / len(densities)
        lower_limit = 0.7 * mean_density
        upper_limit = 1.3 * mean_density
        self.production_run(self.topology_file, self.output_dir)
        if lower_limit <= x <= upper_limit:
            print(f"The given value x = {x} is within 10% accuracy of the mean density y = {mean_density:.2f}.")
            return True
        else:
            print(f"The given value x = {x} is not within 10% accuracy of the mean density y = {mean_density:.2f}.")
            return False


    # def denProd(self,density):
    #     x = density  # Replace with your given density value
    #     is_within_accuracy = self.check_density_accuracy(x, self.densities)
    #
    #     if is_within_accuracy:
    #         self.production_run(self.topology_file, self.output_dir)
    #     else:
    #         print("The density does not satisfy the 10% accuracy condition. The production run will not continue.")

    def correct(self):
        self.trjconv = gromacs.tools.Trjconv(s=self.tpr_file, f=self.trr_file, o=self.xtc_file, pbc="mol", ur="compact")
        self.trjconv.run(input="System")

    # Analysis

    # MSD

    def index_file(self):
        subprocess.run(["gmx_mpi", "make_ndx", "-f", self.gro_file, "-o", f"{self.output_dir}/index.ndx"], input=b"q\n")
    def msd_mkaer(self):
        subprocess.run(
            ["gmx_mpi", "msd", "-f", self.xtc_file, "-s", self.tpr_file, "-n", f"{self.output_dir}/index.ndx", "-o", f"{self.output_dir}/msd.xvg"],
            input=b"0\n")





    def extract_residues_from_itp(self):
        residues = []

        with open(self.nmol_itp, "r") as f:
            lines = f.readlines()
            # for line in lines:
            #     if line.startswith(";"):
            #         continue
            #     if line.strip() == "[ molecules ]":
            #         break

            for line in lines[lines.index("[ molecules ]\n") + 1:]:
                if line.strip() == "" or line.startswith(";"):
                    continue
                residues.append(line.split()[0])

        return residues


    # residue_names = extract_residues_from_itp(nmol_itp)
    # print("Residue names:", residue_names)

    # universes = {}

    # for residue_name in residue_names:
    #     universe = mda.Universe(gro_file)
        # residue_universe = universe.select_atoms(f"resname {residue_name}")
        # universes[residue_name] = residue_universe

    # Usage example:
    # residue_name = "ACN"
    # print(f"Number of '{residue_name}' atoms: {len(universes[residue_name].atoms)}")

    # import MDAnalysis as mda
    # import MDAnalysis.analysis.rdf as rdf
    # import matplotlib.pyplot as plt
    # from matplotlib.backends.backend_pdf import PdfPages

    # universe = mda.Universe(gro_file, xtc_file)


    # Function to count the number of atoms in a residue
    def count_atoms_in_residue(self, universe, resname):
        residue = universe.select_atoms(f"resname {resname}")[0].residue
        return len(residue.atoms)


    # Calculate the number of atoms in each residue type


    # def rdf(self, residue_names):
    #     atoms_in_residues = {resname: self.count_atoms_in_residue(self.universe, resname) for resname in residue_names}
    #     with PdfPages(f'{self.output_dir}/RDF_plots.pdf') as pdf:
    #         # Calculate RDF for each residue type as reference
    #         for reference_residue in residue_names:
    #             plt.figure(figsize=(8, 6))
    #             for target_residue in residue_names:
    #                 reference_group = self.universe.select_atoms(f"resname {reference_residue}")
    #                 target_group = self.universe.select_atoms(f"resname {target_residue}")
    #
    #                 # Use exclusion_block based on the number of atoms in the residues
    #                 exclusion_block = (atoms_in_residues[reference_residue],
    #                                    atoms_in_residues[target_residue]) if reference_residue == target_residue else None
    #
    #                 rdf_analysis = rdf.InterRDF(reference_group, target_group, nbins=75, range=(0.0, 15.0),
    #                                             exclusion_block=exclusion_block)
    #                 rdf_analysis.run()
    #
    #                 # Plot RDF
    #                 plt.plot(rdf_analysis.bins, rdf_analysis.rdf, label=f"{reference_residue}-{target_residue}")
    #
    #             plt.xlabel("Distance (angstroms)")
    #             plt.ylabel("RDF")
    #             plt.title(f"RDF for reference residue: {reference_residue}")
    #             plt.legend()
    #
    #             # Display the plot in the Jupyter notebook
    #             # plt.show()
    #
    #             # Save the plot to the PDF file
    #             pdf.savefig()
    #             plt.close()
    #             return rdf_analysis

    # import numpy as np


    # Function to calculate coordination number from RDF
    def calculate_coordination_number(self,rdf_analysis, rcut):
        # Integrate the RDF up to the cutoff radius
        cn = np.trapz(rdf_analysis.rdf[rdf_analysis.bins <= rcut], x=rdf_analysis.bins[rdf_analysis.bins <= rcut])
        return cn

    def rdf(self, residue_names, rcut):
        universe = mda.Universe(self.gro_file, self.xtc_file)
        atoms_in_residues = {resname: self.count_atoms_in_residue(universe, resname) for resname in residue_names}
        with PdfPages(f'{self.output_dir}/RDF_plots2CordinationNum.pdf') as pdf:
            coordination_numbers = {}
            # Calculate RDF for each residue type as reference
            for reference_residue in residue_names:
                plt.figure(figsize=(8, 6))
                coordination_numbers[reference_residue] = {}
                for target_residue in residue_names:
                    reference_group = universe.select_atoms(f"resname {reference_residue}")
                    target_group = universe.select_atoms(f"resname {target_residue}")

                    # Use exclusion_block based on the number of atoms in the residues
                    exclusion_block = (atoms_in_residues[reference_residue],
                                       atoms_in_residues[target_residue]) if reference_residue == target_residue else None

                    rdf_analysis = rdf.InterRDF(reference_group, target_group, nbins=75, range=(0.0, 15.0),
                                                exclusion_block=exclusion_block)
                    rdf_analysis.run()

                    # Calculate coordination number and store it
                    cn = self.calculate_coordination_number(rdf_analysis, rcut)  # Adjust rcut value as needed
                    coordination_numbers[reference_residue][target_residue] = cn

                    # Plot RDF
                    plt.plot(rdf_analysis.bins, rdf_analysis.rdf, label=f"{reference_residue}-{target_residue} (CN={cn:.2f})")

                plt.xlabel("Distance (angstroms)")
                plt.ylabel("RDF")
                plt.title(f"RDF for reference residue: {reference_residue}")
                plt.legend()

                # Display the plot in the Jupyter notebook
                # plt.show()

                # Save the plot to the PDF file
                pdf.savefig()
                plt.close()
                return coordination_numbers

    # Print coordination numbers
    # for ref_residue, target_residues in coordination_numbers.items():
    #     for target_residue, cn in target_residues.items():
    #         print(f"Coordination number of {ref_residue} to {target_residue}: {cn:.2f}")



    #Create an empty list to store the data
    def cordination_number(self, coordination_numbers):
        coordination_numbers_list = []

        for ref_residue, target_residues in coordination_numbers.items():
            for target_residue, cn in target_residues.items():
                # Append the data to the list instead of printing
                coordination_numbers_list.append([ref_residue, target_residue, cn])

        # Convert the list to a DataFrame
        coordination_numbers_df = pd.DataFrame(coordination_numbers_list,
                                               columns=['Reference_Residue', 'Target_Residue', 'Coordination_Number'])

        # Print the DataFrame
        print(coordination_numbers_df)

    # def available_threads(self):
    #     pass








