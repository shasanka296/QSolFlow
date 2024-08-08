"""
Parse and push data from NLP parsing and submit calculations for extracted molecules
"""
import argparse
import random
from pathlib import Path
from fireworks import LaunchPad
from monty.serialization import dumpfn, loadfn
from d3tales_fw.workflows.wf_writer import *

BASE_DIR = Path(__file__).resolve().parent.parent

parser = argparse.ArgumentParser(description='Launch MD calculations')
parser.add_argument('-f', '--filename', type=str, help='filepath for a JSON nlp data file', default="")
parser.add_argument('-p', '--priority', type=int, help='jobs priority', default=5)
args = parser.parse_args()


def populate_md_wf(path=None,**kwargs):
    lpad_file = os.path.join(BASE_DIR.parent, 'launch', 'md_launchpad.yaml')
    wf = d3tales_md_wf(param_file=path,**kwargs)
    info = LaunchPad().from_file(lpad_file).add_wf(wf)
    fw_id = list(info.values())[0]
    return fw_id


if not os.path.isfile(args.filename):
    systems=["methane_water",'ethane-water']
    key_dic={}
    for i in systems:
       key_dic[i]=random.randint(1,30000000)
    md_kwargs = {
        "smiles_list": ["CO","CO" ,"CCO","CCO"],
        "con_list": [0,0,0,0], "WF_name1":"Test_run1","WF_name2":"Test_run2",  "name_list":["wa1ter","wa2ter","methane1","ethane2"], "cons":0, "type_list":["Solvent","Solvent","Solute1","Solute1"], "dir":"/mnt/gpfs2_4m/scratch/sla296/test_run/output_of_runs", "solvent_name1":["wa1ter"], "solute_name1":["methane"],"solvent_name2":["wa2ter"], "solute_name2":["ethane"], "solvent_smiles1":["CO"], "solute_smiles1":["CCO"],"solvent_smiles2":["CO"], "solute_smiles2":["CCO"],"x1":10, "y1":10, "z1":10, "conmatrix1":["0.1"], "den1":"10", "MM1":"5", "x2":10, "y2":10, "z2":10, "conmatrix2":["0.1"], "den2":"10", "MM2":"5","num_systems":"2", "populate_name":"test","key_dic":key_dic}
    #md_kwargs = {
        #"smiles_list": ["O","C"],
        #"con_list": [0,0], "WF_name1":"Test_run1","WF_name2":"Test_run2",  "name_list":["water","methane1"], "cons":0, "type_list":["Solvent","Solute1"], "dir":"/mnt/gpfs2_4m/scratch/sla296/test_run/output_of_runs", "solvent_name1":["water"], "solute_name1":["methane"],"solvent_name2":["wa2ter"], "solute_name2":["ethane"], "x1":10, "y1":10, "z1":10, "conmatrix1":["0.1"], "den1":"10", "MM1":"10", "x2":10, "y2":10, "z2":10, "conmatrix2":["0.1"], "den2":"10", "MM2":"18.012","num_systems":"1", "populate_name":"test","key_dic":key_dic}

    all_ids = {"test_md_fw": populate_md_wf(**md_kwargs)}

else:
    all_md_data = loadfn(args.filename)
    all_ids = {}
    for mol_name, md_kwargs in all_md_data.items():
        all_ids[mol_name] = populate_md_wf(**md_kwargs)

# dumpfn(all_ids, "md_data/md_{}_molIDs.json".format(args.filename.split("/")[-1]), indent=2)


#md_kwargs = {
        #"smiles_list": ["O","O" ,"C","CC"],
        #"con_list": [0,0,0,0], "WF_name1":"Test_run1","WF_name2":"Test_run2",  "name_list":["wa1ter","wa2ter","methane1","ethane2"], "cons":0, "type_list":["Solvent","Solvent","Solute1","Solute1"], "dir":"/mnt/gpfs2_4m/scratch/sla296/test_run/output_of_runs", "solvent_name1":["wa1ter"], "solute_name1":["methane"],"solvent_name2":["wa2ter"], "solute_name2":["ethane"], "x1":10, "y1":10, "z1":10, "conmatrix1":["0.1"], "den1":"55.508", "MM1":"18.012", "x2":10, "y2":10, "z2":10, "conmatrix2":["0.1"], "den2":"55.508", "MM2":"18.012","num_systems":"2", "populate_name":"test","key_dic":key_dic}
