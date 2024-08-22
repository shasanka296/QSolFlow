import os
import subprocess


class reorg:
    """
    A class for reorganizing the itp files in preparation for simulations.

    Args:
    - molecule (str): The name of the molecule.
    - dir (str): The directory path where reorganized files will be stored.
    - key (str): A key identifier.

    Attributes:
    - molecule (str): The name of the molecule.
    - dir (str): The directory path.
    - moleclue (str): A processed version of the molecule name.

    Methods:
    - __init__(self, molecule, dir, key): Initializes the reorg object, creates necessary files, and reorganizes the ITP files.
    """

    def __init__(self, molecule, direc, key):

        self.moleclue = molecule
        self.dir = direc
        self.atomtype_is_present = True

        subprocess.run(f"touch {os.path.join(self.dir, f'{self.moleclue}.itp')}", shell=True)
        if not os.path.isfile(os.path.join(self.dir, f'{self.moleclue}_atomtype.itp')):
            print("atomtype not found")
            self.atomtype_is_present = False
            subprocess.run(f'touch {os.path.join(self.dir, f"{self.moleclue}_atomtype.itp")}', shell=True)
        with open(f"{os.path.join(self.dir, f'{self.moleclue}f.itp')}", 'r') as org:
            lines = org.readlines()
            self.atomtype = 10
            self.atomtype_lastline = 0
            self.orginal = []
            for iteams in lines:
                self.orginal.append(iteams)
            self.orginal[6] = ";[ defaults ]" + "\n"
            self.orginal[7] = "; nbfunc        comb-rule       gen-pairs       fudgeLJ fudgeQQ" + "\n"
            self.orginal[8] = ";    1               3              yes            0.5     0.5" + "\n"
            for iteams in self.orginal:
                if iteams.strip() == "[ moleculetype ]":
                    self.atomtype_lastline = self.orginal.index(iteams) - 2
            self.write_itp()
            if not self.atomtype_is_present:
                self.write_type()

    def write_itp(self):
        with open(f'{os.path.join(self.dir, f"{self.moleclue}.itp")}', 'a') as itp:
            for i in range(self.atomtype):
                itp.writelines(self.orginal[i])
            for j in range(self.atomtype_lastline + 2, len(self.orginal)):
                itp.writelines(self.orginal[j])

    def write_type(self):
        with open(f'{os.path.join(self.dir, f"{self.moleclue}_atomtype.itp")}', 'a') as a_type:
            a_type.writelines('\n')
            for k in range(self.atomtype_lastline - self.atomtype + 1):
                a_type.writelines(self.orginal[self.atomtype + k])
