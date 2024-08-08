import rdkit.Chem as Chem
from rdkit.Chem import AllChem, Descriptors
import rdkit as rd
import subprocess
import os
from d3tales_fw.workflows.envwf import CONDAPATH


def run_on_terminla(name, q, smiles, direc=None, mul=None):
    charge_mat=[]
    command = (
    f"source {CONDAPATH} && "
    f"conda activate vlxenv && "
    f"python -c \"import d3tales_fw.Fast.new_dft as a; a.get_charge('{name}', {q}, '{smiles}', direc= '{direc}', mul = {mul})\"")

    subprocess.run(command, shell=True,check=True)
    with open (f"{os.path.join(direc,'charge.txt')}", 'r') as file:
        lines= file.readlines()
        for a in lines:
            charge_mat.append(float(a))
    subprocess.run([f"rm {os.path.join(direc,'charge.txt')}"], shell=True,check=True)
    return charge_mat

def get_charge(name, ch, smiles, direc=None, mul=None):
    q= ch if ch else 0
    try:
        import veloxchem as vlx
    except:
        raise Exception("cannot load the module")
    path_to_mol = f"{name}.pdb"
    path_to_xyz = f"{name}.xyz"
    tr_mol = rd.Chem.MolFromSmiles(smiles)  # I found that if I just use the xyz file RDKIT
    # cannot count num of unpaired e, so I am defining the molecule using smiles first
    b = rd.Chem.AddHs(tr_mol)  # adding H so i can later turn pdb
    rd.Chem.AllChem.EmbedMolecule(b)  # Generating 3D rep of molecule
    if direc:
        path_to_mol = os.path.join(direc,name, f"{name}.pdb")
        path_to_xyz = os.path.join(direc,name, f"{name}.xyz")
    print(path_to_xyz)
    molecule = vlx.Molecule.read_xyz_file(path_to_xyz)  # Loading Molecule into veloxchem
    basis = vlx.MolecularBasis.read(molecule, "cc-pVDZ")  # setting up for calc
    unpaired = Descriptors.NumRadicalElectrons(tr_mol)  # finding num of unpaired e for multiplicity calc
    multiplicity=mul
    if mul==None:
        multiplicity = unpaired + 1 + int(q)%2 # 2S+1+charge mod 2
    print(f"the spin multi: {multiplicity}")
    vlx.Molecule.set_multiplicity(molecule, multiplicity)  # setting multi
    vlx.Molecule.set_charge(molecule, q)  # Setting Charge
    scf_drv = vlx.ScfUnrestrictedDriver()  # Unrestricted for sate>singlet
    if multiplicity == 1:
        scf_drv = vlx.ScfRestrictedDriver()  # for singlet
    scf_drv.ostream.mute()  # shutting off the output stream for the calc

    scf_drv.xcfun = "B3LYP"
    # avaible functioals: ['SLATER', 'SLDA', 'B88X', 'BLYP', 'B3LYP', 'BHANDH', 'BHANDHLYP', 'PBE', 'PBE0', 'REVPBE',
    # 'BP86', 'PW91', 'MPW1K', 'OLYP', 'O3LYP', 'X3LYP', 'B97', 'B97-1', 'B97-2', 'B97-3', 'TPSS', 'TPSSH',
    # 'REVTPSS', 'PKZB', 'SCAN', 'RSCAN', 'R2SCAN', 'M05', 'M05-2X', 'M06', 'M06-2X', 'M06-HF', 'M06-L', 'M11-L',
    # 'MPW1B95', 'MPWB1K', 'PW6B95', 'PWB6K']

    scf_drv.grid_level = 4  # grid spacing for numerical calc, i dont know how to pick this but this seems to be the
    # defult
    scf_drv.conv_thresh = 1.0e-6  # setting criteria of deltaE for convergence

    try:
        scf_results = scf_drv.compute(molecule, basis)
        esp_drv = vlx.RespChargesDriver()  #restraied ESP calc driver
        esp_charges = esp_drv.compute(molecule, basis, scf_results, "esp")  # calc
    except AssertionError:
        print(f" Error in DFT clac, please recheck the inputs and try again or manually create the parameters")
        return None
    with open (f"{os.path.join(direc,'charge.txt')}", 'a') as file:
        print(esp_charges.tolist())
        for i in esp_charges.tolist():
            file.write(f"{str(i)}\n")
    return None

