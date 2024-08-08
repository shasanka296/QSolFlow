import argparse
from pathlib import Path
from fireworks import LaunchPad
from d3tales_fw.workflows.wf_writer import *
from d3tales_fw.workflows.ParamSet import GausParamSet

BASE_DIR = Path(__file__).resolve().parent.parent


def populate_d3tales_lpad(
    filename=None,
    ids_list=None,  # norm: None (must be None for json to be read)
    smiles_list=None,  # norm: None  (must be None for json to be read)
    required_name_string='',  # norm: ''

    workflow_function=d3tales_wf,
    priority=4,

    public=False,
    default_origin='Zinc',
    submit=True,  # norm: True
    check_if_already_run=True,  # norm: True
    solvent='acetonitrile',  # norm: acetonitrile  ( DiMethylSulfoxide )
    param_tag='gaus_',  # norm: gaus_  ( huckaba_ dihedral_ gaus_singlets_)
    hf_mol_opt=False,  # norm: False
    restricted=True,  # norm: True
    use_iop=True,  # norm: True
    wtune=True,  # norm: True
    run_nto=False,  # norm: False

    email=None,  # "rdu230@uky.edu"
    username=None,  # "Rebekah Duke"
    wf_tag=""  # ""
):
    """
    Function to initiate a series of D3TaLES workflows for a given set of molecules.
    """
    gaussian_file_name = 'gaussian_ur' if not restricted else 'gaussian_noTune' if not use_iop else 'gaussian_hf' if hf_mol_opt else 'gaussian'
    lpad_file = os.path.join(BASE_DIR.parent, 'launch', 'md_launchpad.yaml')
    param_file = os.path.join(BASE_DIR, 'parameters', param_tag + 'parameter_file.json')

    smiles_dict = {}
    if filename and os.path.isfile(filename):
        with open(filename, 'r') as f:
            smiles_dict = json.load(f)

    origin = smiles_dict.pop("origin") if smiles_dict.get("origin") else default_origin
    paramset = GausParamSet().from_json(param_file)
    mol_list = ids_list or smiles_list or list(smiles_dict.values())
    names_dict = {smiles: name for name, smiles in smiles_dict.items()}

    fw_id_list = []
    kwarg_dict = dict(paramset=paramset, public=public, hf_mol_opt=hf_mol_opt, submit=submit, use_iop=use_iop,
                      check_if_already_run=check_if_already_run, solvent=solvent, priority=priority,
                      gaussian_file_name=gaussian_file_name, email=email, username=username, restricted=restricted,
                      wtune=wtune, run_nto=run_nto, origin_group=origin, wf_tag=wf_tag)
    for mol in mol_list:
        if ids_list:
            kwarg_dict.update(identifier=mol)
        else:
            name = names_dict.get(mol)
            if required_name_string not in name:
                continue
            kwarg_dict.update(smiles=mol, mol_name=name)
        wf = workflow_function(**kwarg_dict)
        info = LaunchPad().from_file(lpad_file).add_wf(wf)
        fw_id = list(info.values())[0]
        fw_id_list.append(fw_id)

    print("Added {} workflows to the launchpad".format(len(mol_list)))
    print("FW ID List: ", ''.join(str(fw_id_list).split(' ')))
    return fw_id_list


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Launch D3TaLES FireWorks workflow.')
    parser.add_argument('filename', metavar='filename', type=str, help='filepath for a JSON molecule file', default=False)
    parser.add_argument('-p', '--priority', type=int, help='jobs priority', default=2)
    parser.add_argument('-l', '--id_list', action='store_true', help='denotes that the filename argument is actually a'
                                                                     'comma seperated list of molecule IDs to submit')
    args = parser.parse_args()

    if args.id_list:
        populate_d3tales_lpad(ids_list=args.filename.split(","), priority=args.priority)
    else:
        populate_d3tales_lpad(filename=args.filename)
