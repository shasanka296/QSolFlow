import re
import json
import argparse
import functools
import pandas as pd
from d3tales_api.D3database.d3database import FrontDB

mol_i_props = ["smiles", "source_group", "groundState_charge", "number_of_atoms", "molecular_weight", "sa_score"]
mol_props = ["hole_reorganization_energy", "electron_reorganization_energy",  "relaxation_groundState_cation1", "relaxation_cation1_groundState", "relaxation_groundState_anion1", "relaxation_anion1_groundState", "vertical_ionization_energy", "vertical_electron_affinity", "adiabatic_ionization_energy", "adiabatic_electron_affinity", "oxidation_potential", "reduction_potential", "rmsd_groundState_cation1", "rmsd_cation1_cation2", "rmsd_groundState_anion1", "rmsd_anion1_anion2", "omega"]
species_props = ["charge", "globular_volume", "radical_stability_score", "homo_lumo_gap", "dipole_moment", "solvation_energy", "homo", "lumo"]

MOL_INFO = ",".join([f"mol_info.{p}" for p in mol_i_props])
MOL_CHAR = ",".join([f"mol_characterization.{p}.0.value" for p in mol_props])
SPECIES_CHAR = ",".join([f"species_characterization.groundState.{p}.0.value" for p in species_props]
                        + [f"species_characterization.cation1.{p}.0.value" for p in species_props]
                        + [f"species_characterization.cation2.{p}.0.value" for p in species_props]
                        + [f"species_characterization.anion1.{p}.0.value" for p in species_props]
                        + [f"species_characterization.anion2.{p}.0.value" for p in species_props]
                        )
ALL_PROPERTIES = ",".join(["public", MOL_INFO, MOL_CHAR, SPECIES_CHAR])


def rgetkeys(_dict, keys, **kwargs):
    def _getkey(_dict, key):
        _dict = _dict or {}
        if isinstance(_dict, dict):
            return _dict.get(key, **kwargs)
        if isinstance(_dict, list) and key.isdigit():
            return _dict[int(key)]

    return functools.reduce(_getkey, [_dict] + keys.split('.'))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Generate a CSV file containing given properties (default: smiles) for all D3TaLES molecules.')
    parser.add_argument('-f', '--filename', type=str, help='filepath for a CSV file to which to save the data',
                        default='structure_data/d3tales.csv')
    parser.add_argument('-p', '--collect_properties', type=str, help='properties to include in CSV',
                        default=ALL_PROPERTIES)
    parser.add_argument('-ids', '--ids_filename', type=str,
                        help='filepath for JSON file containing ids (as keys) to search')
    args = parser.parse_args()

    # Get property name and generate query
    collect_properties = args.collect_properties.split(',')
    prop_paths = {p.split(".0.")[0].replace('mol_info.', "").replace('species_characterization.', "").replace(
        'mol_characterization.', ""): p for p in collect_properties}
    projection = {p.split(".0.")[0].strip('.'): 1 for p in collect_properties}
    print(','.join(prop_paths.keys()))

    # Database search
    frontend_db = FrontDB(collection_name="base")
    if args.ids_filename:
        with open(args.ids_filename) as f:
            ids_data = json.load(f)
        ids = list(ids_data.keys())
        cursor = frontend_db.coll.find({'_id': {"$in": ids}}, {projection: 1})
    else:
        cursor = frontend_db.coll.find({}, projection)

    # Data cleaning with pandas
    master_data = pd.DataFrame.from_records(cursor)
    master_data.set_index('_id', inplace=True)
    print(master_data.head())
    for prop_name, prop_path in prop_paths.items():
        master_data[prop_name] = master_data.apply(lambda x: rgetkeys(x.to_dict(), prop_path), axis=1)
    final_data = master_data[prop_paths.keys()]
    final_data.to_csv(args.filename)
    print("Sucess! D3TaLES molecules and {} property were saved to {}".format(args.collect_properties, args.filename))
