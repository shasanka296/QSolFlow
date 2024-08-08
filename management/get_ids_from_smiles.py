import os
import json
import argparse
from d3tales_api.D3database.d3database import FrontDB
from rdkit.Chem import MolFromSmiles, MolToSmiles

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Get D3TaLES IDs from SMILES JSON file.')
    parser.add_argument('-f', '--filename', type=str, help='filepath for a JSON file from which to get SMILES', default='')
    parser.add_argument('-i', '--ids_list', type=str, help='list of ids from which to get smiles', default='')
    parser.add_argument('-g', '--generate_id', help='generate D3TaLES ID if none exists for given SMILES', action='store_true')
    parser.add_argument('-s', '--return_smiles', help='return SMILES as value instead of name', action='store_true')
    args = parser.parse_args()

    if os.path.isfile(args.filename):
        with open(args.filename, 'r') as f:
            smiles_dict = json.load(f)
        origin = smiles_dict.pop("origin") if smiles_dict.get("origin") else 'Risko'

        id_data = {}
        for name, smiles in smiles_dict.items():
            value = smiles if args.return_smiles else name
            try:
                rdkmol = MolFromSmiles(smiles)
                clean_smiles = MolToSmiles(rdkmol)
            except:
                continue
            if args.generate_id:
                _id = FrontDB(smiles=clean_smiles, group=origin).generate_id()
                id_data[_id] = value
            else:
                _id = FrontDB(smiles=clean_smiles, group=origin).check_if_in_db()
                if _id:
                    id_data[_id] = value
        with open("structure_data/ids_from_smiles.json", 'w') as f:
            json.dump(id_data, f, indent=2)
        print(id_data)
        print(','.join(id_data.keys()))

    if args.ids_list:
        return_dict = {}
        for i in args.ids_list.split(','):
            query = FrontDB().make_query({"_id": i}, {"mol_info.smiles": 1})
            return_dict[i] = query.mol_info[0].get("smiles")
        with open("structure_data/ids_smiles.json", 'w') as f:
            json.dump(return_dict, f, indent=2)





