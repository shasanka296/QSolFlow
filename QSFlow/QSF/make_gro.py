import subprocess
import os
import QSFlow.QSF.reorg as g
import QSFlow.QSF.topology as top
from QSFlow.workflows.envwf import MDP_Location


class Gro:
    """
    A class to manage Gro files out of the packmol output and clean up the file system to feed to the simulation code.
        Methods:
            - initial(): Moves itps and initial files to run the simulation code, and cleans up the file system for
              simulation.

    """
    def __init__(self, solvent, solute, solvent2, di, x, y, z, key, initla_system, path_to_file):
        """

        :param solvent: The full name of the solvent.
        :param solute: The full name of the solute.
        :param solvent2: The full name of the second solvent. For future use.
        :param di: Directory where QSF outputs are located.
        :param x: X-length of the topology.
        :param y: Y-length of the topology.
        :param z: Z-length of the topology.
        :param key: Key to the current system.
        :param initla_system: Boolean to specify if initial params are given.
        :param path_to_file: Path to inital files.


        """
        self.dir = di
        self.x = x
        self.y = y
        self.z = z
        self.fullname_solvent = solvent
        self.solvent = solvent[:3]
        self.key = key
        self.command12 = f'cp -r {MDP_Location} {os.path.join(self.dir, f"InputGrofiles{self.key}")}'
        self.solute = solute
        self.solvent2 = solvent2
        self.initla_system = initla_system
        self.path_to_file = path_to_file
        if self.initla_system:
            self.initial()
        else:

            print(f"in the gro solutes: {self.solute}")

            command_main = f'gmx_mpi editconf -f {os.path.join(self.dir, f"Packmol{key}", "Mixture.pdb")} -box {self.x} {self.y} {self.z} -o {os.path.join(self.dir, f"Packmol{key}", "solvated.gro")}'
            command_alt = f'gmx editconf -f {os.path.join(self.dir, f"Packmol{key}", "Mixture.pdb")} -box {self.x} {self.y} {self.z} -o {os.path.join(self.dir, f"Packmol{key}", "solvated.gro")}'
            try:
                subprocess.run(command_main, shell=True, check=True)
            except subprocess.CalledProcessError:
                print("Looks like there is no MPI")
                subprocess.run(command_alt, shell=True, check=True)

            print("gro file made")

            command4 = f'cp {os.path.join(self.dir, f"{solvent[:3]}_Solvent", f"{solvent[:3]}_Solvent.pdb")} {os.path.join(self.dir, f"InputGrofiles{key}")}'

            command5 = f'cp {os.path.join(self.dir, f"Packmol{key}", "solvated.gro")} {os.path.join(self.dir, f"InputGrofiles{key}")}'

            command6 = f'rm -r {os.path.join(self.dir, f"Packmol{key}")}'

            command8 = f'cp {os.path.join(self.dir, f"{solvent2}_Solvent2", f"{solvent2}_Solvent2.pdb")} {os.path.join(self.dir, "InputGrofiles")} && mv {os.path.join(self.dir, f"{solvent2}_Solvent2", f"{solvent2}_Solvent2.gmx.itp")} {os.path.join(self.dir, f"InputGrofiles{key}")}'

            command10 = f'rm -r {os.path.join(self.dir, f"{solvent2}_Solvent2")}'

            print("cleaning up")

            subprocess.run(command4, shell=True, check=True)

            subprocess.run(command5, shell=True, check=True)

            subprocess.run(command6, shell=True, check=True)

            if len(self.solvent2) != 0:
                subprocess.run(command8, shell=True, check=True)

                subprocess.run(command10, shell=True, check=True)

            for i in range(len(self.solute)):
                com1 = f"cp {os.path.join(self.dir, f'{self.solute[i].strip()[:3]}_Solute{1}', f'{self.solute[i].strip()[:3]}_Solute{1}.pdb')} {os.path.join(self.dir, f'InputGrofiles{key}')}"

                subprocess.run(com1, shell=True, check=True)

            g.reorg(self.solvent + "_Solvent", os.path.join(self.dir, f"InputGrofiles{key}"), key)

            for j in self.solute:
                print(f'ran {j}')

                g.reorg(j.strip()[:3] + f"_Solute1", os.path.join(self.dir, f"InputGrofiles{key}"), key)

            subprocess.run(self.command12, shell=True, check=True)

        print("Starting the simulation")

    def initial(self):

        top.toopol(currentdir=self.dir, key=self.key, inital=True)

        try:

            self.commandMDP = f"cp -r {os.path.join(self.path_to_file, 'MDP')} {os.path.join(self.dir, f'InputGrofiles{self.key}')}"

            subprocess.run(self.commandMDP, shell=True, check=True)

        except subprocess.CalledProcessError:

            subprocess.run(self.command12, shell=True, check=True)

        move_the_initial_system = f"cp {os.path.join(self.path_to_file, f'{self.fullname_solvent}_{self.solute[0]}.gro')} {os.path.join(self.dir, f'InputGrofiles{self.key}', 'solvated.gro')}"

        move_topol = f"cp {os.path.join(self.path_to_file, 'topol.top')} {os.path.join(self.dir, f'InputGrofiles{self.key}')}"

        move_nmol = f"cp {os.path.join(self.path_to_file, 'nmol.itp')} {os.path.join(self.dir, f'InputGrofiles{self.key}')}"

        move_solvent_itp = f"cp {os.path.join(self.path_to_file, f'{self.fullname_solvent}.itp')} {os.path.join(self.dir, f'InputGrofiles{self.key}')}"

        move_solvent_atomtype = f"cp {os.path.join(self.path_to_file, f'{self.fullname_solvent}_atomtypes.itp')} {os.path.join(self.dir, f'InputGrofiles{self.key}')}"

        subprocess.run(move_the_initial_system, shell=True, check=True)

        subprocess.run(move_topol, shell=True, check=True)

        subprocess.run(move_nmol, shell=True, check=True)

        subprocess.run(move_solvent_itp, shell=True, check=True)

        try:

            subprocess.check_call(move_solvent_atomtype, shell=True)

            subprocess.run(move_solvent_atomtype, shell=True, check=True)

        except subprocess.CalledProcessError:

            print("Looks like you did not provide the atomtypes for your solvent. The system continued without it.")

        for i in self.solute:

            move_solute_itp = f"cp {os.path.join(self.path_to_file, f'{i}.itp')} {os.path.join(self.dir, f'InputGrofiles{self.key}')}"

            move_solute_atomtype = f"cp {os.path.join(self.path_to_file, f'{i}_atomtypes.itp')} {os.path.join(self.dir, f'InputGrofiles{self.key}')}"

            subprocess.run(move_solute_itp, shell=True, check=True)

            try:

                subprocess.check_call(move_solute_atomtype, shell=True)

                subprocess.run(move_solute_atomtype, shell=True, check=True)

            except subprocess.CalledProcessError:

                print(f"Looks like you did not provide the atomtypes for {i}. The system continued without it.")
