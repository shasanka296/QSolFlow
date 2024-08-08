import subprocess
import os
import d3tales_fw.Fast.reorg as g
import d3tales_fw.Fast.topology as top
from d3tales_fw.workflows.envwf import MDP_Location,meta_dir,CONDAPATH

class gro:
    def __init__(self, solvent, solute, solvent2, di, x, y, z, key, initla_system,path_to_file):
        self.dir = di
        self.x = x
        self.y = y
        self.z = z
        self.fullname_solvent=solvent
        self.solvent = solvent[:3]
        self.key=key
        self.command12 = f'cp -r {MDP_Location} {os.path.join(self.dir,f"InputGrofiles{self.key}")}'
        self.solute = solute
        self.solvent2 = solvent2
        self.initla_system = initla_system
        self.path_to_file = path_to_file
        if self.initla_system:
            self.initial()
        else:

            print(f"in the gro solutes: {self.solute}")

            command = f'gmx_mpi editconf -f {os.path.join(self.dir, f"Packmol{key}", "Mixture.pdb")} -box {self.x} {self.y} {self.z} -o {os.path.join(self.dir, f"Packmol{key}", "solvated.gro")}'

            subprocess.run(command, shell=True, check=True)

            print("gro file made")

            command4 = f'mv {os.path.join(self.dir, f"{solvent[:3]}_Solvent", f"{solvent[:3]}_Solvent.pdb")} {os.path.join(self.dir, f"InputGrofiles{key}")}'

            command5 = f'mv {os.path.join(self.dir, f"Packmol{key}", "solvated.gro")} {os.path.join(self.dir, f"InputGrofiles{key}")}'

            command6 = f'rm -r {os.path.join(self.dir, f"{solvent[:3]}_Solvent")} && rm -r {os.path.join(self.dir, f"Packmol{key}")}'

            command8 = f'mv {os.path.join(self.dir, f"{solvent2}_Solvent2", f"{solvent2}_Solvent2.pdb")} {os.path.join(self.dir, "InputGrofiles")} && mv {os.path.join(self.dir, f"{solvent2}_Solvent2", f"{solvent2}_Solvent2.gmx.itp")} {os.path.join(self.dir, f"InputGrofiles{key}")}'

            command10 = f'rm -r {os.path.join(self.dir, f"{solvent2}_Solvent2")}'

            print("cleaning up")

            subprocess.run(command4, shell=True, check=True)

            subprocess.run(command5, shell=True, check=True)

            subprocess.run(command6, shell=True, check=True)

            if len(self.solvent2) != 0:
                subprocess.run(command8, shell=True, check=True)

                subprocess.run(command10, shell=True, check=True)

            for i in range(len(self.solute)):
                com1 = f"mv {os.path.join(self.dir, f'{self.solute[i].strip()[:3]}_Solute{1}', f'{self.solute[i].strip()[:3]}_Solute{1}.pdb')} {os.path.join(self.dir, f'InputGrofiles{key}')}"

                com2 = f'rm -r {os.path.join(self.dir, f"{self.solute[i].strip()[:3]}_Solute{1}")}'

                subprocess.run(com1, shell=True, check=True)

                subprocess.run(com2, shell=True, check=True)

            g.reorg(self.solvent + "_Solvent", os.path.join(self.dir, f"InputGrofiles{key}"), key)

            if len(self.solvent2) != 0:
                g.reorg(self.solvent2 + "_Solvent2", os.path.join(self.dir, f"InputGrofiles{key}"))

            for j in range(len(self.solute)):
                print(f'ran solute{j + 1}')

                g.reorg(self.solute[j].strip()[:3] + f"_Solute{j + 1}", os.path.join(self.dir, f"InputGrofiles{key}"), key)

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
