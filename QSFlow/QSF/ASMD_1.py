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
from QSFlow.workflows.envwf import meta_dir


class ASMD:
    """A class that uses the GROMACS wrapper along with subprocess to create and run simulation using GROMACS 2020
    software package to find energies and uses pandas, NumPy and matplotlib to store and visualize data, as well as
    to generate  RDF plots. A density check function is built into the class though, it is disabled due to lack of
    reliable data about solution densities for arbitrary mixtures. The order of function class for the standard
    simulation is as follows: count_warnings, EnergyMin, NVT, NPT, calculate_density, check_density_accuracy (NOTE:
    CURRENTLY THE FUNCTION IS CALLED FOR THE SAKE OF CONSISTENCY, IT DOES NOT ACTUALLY VERIFY THE DENSITY), correct,
    index_file, extract_residues_from_itp, RDF, coordination_number.

    """

    def __init__(self, key):

        self.current_dir = meta_dir
        self.input_dir = os.path.join(self.current_dir, f"InputGrofiles{key}")
        self.output_dir = os.path.join(self.current_dir, f"Output{key}")
        self.mdp_dir = os.path.join(self.current_dir, f"InputGrofiles{key}", "MDP")
        self.topology_file = os.path.join(self.input_dir, "topol.top")
        self.initial_coordinates = os.path.join(self.input_dir, "solvated.gro")
        self.xtc_file = os.path.join(self.output_dir, f"production.xtc")
        self.tpr_file = os.path.join(self.output_dir, f"production.tpr")
        self.eqtpr_file = os.path.join(self.output_dir, f"equilibration.tpr")
        self.trr_file = os.path.join(self.output_dir, f"production.trr")
        self.nmol_itp = os.path.join(self.input_dir, "nmol.itp")
        self.gro_file = os.path.join(self.output_dir, f"production.gro")
        self.edr_file = os.path.join(self.output_dir, f"equilibration.edr")
        self.MPI = True
        try:
            subprocess.run(['gmx_mpi'], shell=True, check=True)
        except subprocess.CalledProcessError:
            self.MPI = False

        self.GMX_prefix = "gmx_mpi" if self.MPI else "gmx"

        subprocess.run(["cp", "-r", f"{self.mdp_dir}/.", self.output_dir])
        os.makedirs(self.output_dir, exist_ok=True)
        self.max_warn = 10
        self.threshold = 100

    def count_warnings(self, filename):
        # Read the simulation output log file and count the number of warnings
        if os.path.isfile(self.output_dir + "/em.log"):
            print("log is made")
            with open(filename, 'r') as file:
                log_data = file.read()

            warning_pattern = r'WARNING'
            warning_count = len(re.findall(warning_pattern, log_data))
        else:
            warning_count = 100

        return warning_count

    @staticmethod
    def get_available_cpu_threads():
        cpu_info = os.popen('lscpu').read()
        threads_pattern = r'Thread\(s\) per core:\s+(\d+)'
        result = re.search(threads_pattern, cpu_info)
        if result:
            threads_per_core = int(result.group(1))
            available_threads = os.cpu_count() * threads_per_core
            return available_threads
        else:
            return None

    @staticmethod
    def get_available_gpus():
        # Get the number of available GPUs
        gpu_info = os.popen('nvidia-smi --list-gpus').read()
        gpu_count = len(re.findall(r'GPU\s\d:', gpu_info))
        return gpu_count

    def run_gromacs_simulation(self, command):
        # Execute the GROMACS simulation command
        try:
            subprocess.run(f'{command} -maxwarn {str(self.max_warn)}', shell=True,
                           check=True, capture_output=True, text=True)
        except subprocess.CalledProcessError as E:
            print("too many errors,updating maxwarn")
            f = E.stderr.split("\n")
            for i in f:
                if "There were" in i:
                    error_line = i.split()
                    self.max_warn = int(error_line[2]) + 10
                    break
            subprocess.run(f'{command} -maxwarn {str(self.max_warn)}', shell=True,
                           check=True, capture_output=True, text=True)

        return "updated maxwarn"

    def EnergyMin(self):
        command = (
            f"{self.GMX_prefix} grompp "
            f"-f {os.path.join(self.mdp_dir, 'em.mdp')} "
            f"-c {self.initial_coordinates} "
            f"-p {self.topology_file} "
            f"-o {os.path.join(self.output_dir, 'em.tpr')}"
        )

        print(command)
        self.run_gromacs_simulation(command)

        available_threads = self.get_available_cpu_threads()
        # available_gpus = self.get_available_gpus()

        command = [f'{self.GMX_prefix}', "mdrun", "-deffnm", "em", "-v"]
        if available_threads is not None:
            command += ["-ntomp", str(available_threads)]
        subprocess.run(command, cwd=self.output_dir)

    def NVT(self):
        command = (
            f"{self.GMX_prefix} grompp "
            f"-f {os.path.join(self.mdp_dir, 'nvt.mdp')} "
            f"-c {os.path.join(self.output_dir, 'em.gro')} "
            f"-p {self.topology_file} "
            f"-o {os.path.join(self.output_dir, 'nvt.tpr')}"
        )

        self.run_gromacs_simulation(command)

        command = [f'{self.GMX_prefix}', "mdrun", "-deffnm", "nvt", "-v"]
        if self.get_available_cpu_threads() is not None:
            command += ["-ntomp", str(self.get_available_cpu_threads())]
        subprocess.run(command, cwd=self.output_dir)

    def NPT(self):
        command = (
            f"{self.GMX_prefix} grompp "
            f"-f {os.path.join(self.mdp_dir, 'equilibration.mdp')} "
            f"-c {os.path.join(self.output_dir, 'nvt.gro')} "
            f"-p {self.topology_file} "
            f"-o {os.path.join(self.output_dir, 'equilibration.tpr')}"
        )

        self.run_gromacs_simulation(command)
        command = [f'{self.GMX_prefix}', "mdrun", "-deffnm", "equilibration", "-v"]
        # if available_gpus > 0:
        # command = ["mpirun", "-np", str(available_gpus), f'{self.GMX_prefix}', "mdrun", "-deffnm", "equilibration", "-v", "-nb", "gpu"]
        if self.get_available_cpu_threads() is not None:
            command += ["-ntomp", str(self.get_available_cpu_threads())]
        subprocess.run(command, cwd=self.output_dir)

    def calculate_density(self):
        gromacs.config.setup()

        density_xvg_file = f"{self.output_dir}/density.xvg"
        try:
            energy = gromacs.tools.Energy_mpi(s=os.path.abspath(self.eqtpr_file),
                                              f=os.path.abspath(self.edr_file),
                                              o=f'{density_xvg_file}')
        except AttributeError:
            energy = gromacs.tools.Energy(s=os.path.abspath(self.eqtpr_file),
                                          f=os.path.abspath(self.edr_file),
                                          o=f'{density_xvg_file}')

        energy.run(input="Density")

        densities = []
        with open(f"{self.output_dir}/density.xvg", "r") as f:
            for line in f.readlines():
                if line.startswith("#") or line.startswith("@"):
                    continue
                densities.append(float(line.split()[1]))

        return densities

    def extract_density(self):
        densities = []
        with open(f'{os.path.join(self.output_dir, "density.xvg")}', "r") as f:
            for line in f.readlines():
                if line.startswith("#") or line.startswith("@"):
                    continue
                densities.append(float(line.split()[1]))

        return densities

    def utilize(self):
        densities = self.extract_density()
        print("Density values:", densities)
        return densities

    def production_run(self):

        command = (
            f"{self.GMX_prefix} grompp "
            f"-f {os.path.join(self.mdp_dir, 'production.mdp')} "
            f"-c {os.path.join(self.output_dir, 'equilibration.gro')} "
            f"-p {self.topology_file} "
            f"-o {os.path.join(self.output_dir, 'production.tpr')}"
        )

        self.run_gromacs_simulation(command)
        command = [f'{self.GMX_prefix}', "mdrun", "-deffnm", "production", "-v"]

        if self.get_available_cpu_threads() is not None:
            command += ["-ntomp", str(self.get_available_cpu_threads())]
        subprocess.run(command, cwd=self.output_dir)
        print("Simulation completed.")

    def check_density_accuracy(self, x, densities):
        mean_density = sum(densities) / len(densities)
        lower_limit = 0.7 * mean_density
        upper_limit = 1.3 * mean_density
        self.production_run()
        if lower_limit <= x <= upper_limit:
            print(f"The given value x = {x} is within 10% accuracy of the mean density y = {mean_density:.2f}.")
            return True
        else:
            print(f"The given value x = {x} is not within 10% accuracy of the mean density y = {mean_density:.2f}.")
            return False

    def correct(self):
        trjconv = gromacs.tools.Trjconv(s=self.tpr_file, f=self.trr_file, o=self.xtc_file, pbc="mol", ur="compact")
        trjconv.run(input="System")

    def index_file(self):
        subprocess.run([f'{self.GMX_prefix}', "make_ndx", "-f", self.gro_file, "-o", f"{self.output_dir}/index.ndx"],
                       input=b"q\n")

    def msd_maker(self):
        subprocess.run(
            [f'{self.GMX_prefix}', "msd", "-f", self.xtc_file, "-s", self.tpr_file, "-n",
             f"{self.output_dir}/index.ndx", "-o", f"{self.output_dir}/msd.xvg"],
            input=b"0\n")

    def extract_residues_from_itp(self):
        residues = []

        with open(self.nmol_itp, "r") as f:
            lines = f.readlines()

            for line in lines[lines.index("[ molecules ]\n") + 1:]:
                if line.strip() == "" or line.startswith(";"):
                    continue
                residues.append(line.split()[0])

        return residues

    @staticmethod
    def count_atoms_in_residue(universe, resname):
        print(resname)
        residue = universe.select_atoms(f"resname {resname}")[0].residue
        return len(residue.atoms)

    @staticmethod
    def calculate_coordination_number(rdf_analysis, rcut):
        cn = np.trapz(rdf_analysis.rdf[rdf_analysis.bins <= rcut], x=rdf_analysis.bins[rdf_analysis.bins <= rcut])
        return cn

    def rdf(self, residue_names, rcut):
        universe = mda.Universe(self.gro_file, self.xtc_file)
        atoms_in_residues = {resname: self.count_atoms_in_residue(universe, resname) for resname in residue_names}
        print(atoms_in_residues)

        with PdfPages(f'{os.path.join(self.output_dir,"RDF_plots2CoordinationNum.pdf")}') as pdf:
            coordination_numbers = {}
            for reference_residue in residue_names:
                plt.figure(figsize=(8, 6))
                coordination_numbers[reference_residue] = {}
                for target_residue in residue_names:
                    reference_group = universe.select_atoms(f"resname {reference_residue}")
                    target_group = universe.select_atoms(f"resname {target_residue}")

                    exclusion_block = (atoms_in_residues[reference_residue],
                                       atoms_in_residues[
                                           target_residue]) if reference_residue == target_residue else None

                    rdf_analysis = rdf.InterRDF(reference_group, target_group, nbins=75, range=(0.0, 15.0),
                                                exclusion_block=exclusion_block)
                    rdf_analysis.run()

                    cn = self.calculate_coordination_number(rdf_analysis, rcut)
                    coordination_numbers[reference_residue][target_residue] = cn

                    plt.plot(rdf_analysis.bins, rdf_analysis.rdf,
                             label=f"{reference_residue}-{target_residue} (CN={cn:.2f})")

                plt.xlabel("Distance (angstroms)")
                plt.ylabel("RDF")
                plt.title(f"RDF for reference residue: {reference_residue}")
                plt.legend()

                pdf.savefig()
                plt.close()
                return coordination_numbers

    @staticmethod
    def coordination_number(coordination_numbers):
        coordination_numbers_list = []

        for ref_residue, target_residues in coordination_numbers.items():
            for target_residue, cn in target_residues.items():
                coordination_numbers_list.append([ref_residue, target_residue, cn])

        coordination_numbers_df = pd.DataFrame(coordination_numbers_list,
                                               columns=['Reference_Residue', 'Target_Residue', 'Coordination_Number'])

        print(coordination_numbers_df)
