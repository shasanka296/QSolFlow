import os
import subprocess
import time
import rdkit as rd
from d3tales_fw.workflows.envwf import BOSSDIR, CONDAPATH,SINGPATH

class lig:
    def __init__(self, smiles, regular_name, molecule, charge, di, own, own_path):
        self.dir = di
        self.smiles = smiles
        self.mol = molecule
        self.charge = int(charge) or 0
        self.path_to_pdb= os.path.join(self.dir, f"{self.mol}.pdb")
        print(self.path_to_pdb)
        path_to_output=os.path.join(self.dir, self.mol)
        if own == True:
            self.own(smiles, regular_name, molecule, charge, dir, own, own_path)
        else:
            tr_mol = rd.Chem.MolFromSmiles(smiles)
            b = rd.Chem.AddHs(tr_mol)
            rd.Chem.AllChem.EmbedMolecule(b)
            rd.Chem.MolToPDBFile(b,self.path_to_pdb)
            conda_activate = f"source {CONDAPATH} && conda activate ligpg"
            export_bossdir = f" export BOSSdir={BOSSDIR}"
            ligpargen_cmd = f"ligpargen -i {self.path_to_pdb} -n { self.mol} -p {path_to_output} -r {regular_name[:3]} -c {self.charge} -o 0 -cgen CM1A"
            singularity_container = SINGPATH
            obabel_cmd=f"obabel -ipdb {self.path_to_pdb} -oxyz -O {os.path.join(self.dir,  self.mol, f'{self.mol}.xyz')} "
            moving_pdb_command= f"mv {self.path_to_pdb} {os.path.join(self.dir,  self.mol )}"

            cmd = ["singularity", "exec", singularity_container, "bash", "-c",
                   f'{conda_activate} && {export_bossdir} && {ligpargen_cmd} && {obabel_cmd} &&  {moving_pdb_command}']
            try:
                print("in try")
                print(cmd)
                subprocess.Popen(cmd).wait()

            except:
                print(f'Ligpargen was not able to find a parameter, user input files is being used. Please rerun with your own itp and pdb files. This is passed for smiles, regular_name, molecule, charge, dir {(smiles, regular_name, molecule, charge, dir)}')
            self.PDBMAKER(self.mol,self.smiles)

    def PDBMAKER(self, name, smiles):
        if name == "MET":
           name = "MeT"
        with open(
                f"{os.path.join(self.dir,self.mol,f'{self.mol}.pdb')}",'r') as fie, open(
                f"{os.path.join(self.dir,self.mol,'new.pdb')}",
                'a') as new:
            f = fie.readlines()
            line_to_print = []
            for iteams in f:
                line_to_print.append(iteams.strip("\n").split(' '))
            for lines in line_to_print:
                if 'UNL' in lines:
                    lines[lines.index('UNL')] = f"{name[:3]}"
            for lines in line_to_print:
                for line in lines:
                    if line == "":
                        new.write(f'')
                    new.write(f'{line} ')
                new.write("\n")
            subprocess.run([f"rm {os.path.join(self.dir,self.mol,f'{self.mol}.pdb')} && mv {os.path.join(self.dir,self.mol,'new.pdb')} {os.path.join(self.dir,self.mol,f'{self.mol}.pdb')}"], shell=True, check=True)
    def own(self,smiles, regular_name, molecule, charge, dir, own, own_path):
        subprocess.run([f'mkdir {os.path.join(self.dir,self.mol)}'], shell=True)
        path1=os.path.join(own_path, f"{regular_name}.pdb")
        path2=os.path.join(own_path, f"{regular_name}.itp")
        cmd1= f" cp  {path1} {os.path.join(self.dir, self.mol, f'{molecule}.pdb')}"
        cmd2 = f" cp {path2} {os.path.join(self.dir, self.mol, f'{molecule}.gmx.itp')} "

        subprocess.run([cmd1], shell=True)
        subprocess.run([cmd2], shell=True)



