import json
import time
import argparse
from d3tales_api.D3database.d3database import FrontDB, D3Database

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Get unfinished IDs in the D3TaLES FireWorks workflow.')
    parser.add_argument('-m', '--mol_characterization', type=str, help='mol_characterization query for finding mol IDs')
    parser.add_argument('-sp', '--special', action='store_true', help='use special query for finding mol IDs')
    parser.add_argument('-fw_sp', '--fw_special', action='store_true', help='use special Fireworks query for finding mol IDs')
    parser.add_argument('-s', '--state', type=str, help='recover FW IDs of a specific type')
    parser.add_argument('-u', '--unfinsihed', action='store_true', help='recover database IDs for molecules with unfinished workflows')
    parser.add_argument('-w', '--workflows', action='store_true', help='search workflow database instead of fireworks database')
    parser.add_argument('-ni', '--not_in', type=str, help='JSON file containing ID keys that should not be included in final ID list', default='')
    args = parser.parse_args()

    frontend_db = FrontDB(collection_name="base")
    fireworks_db = D3Database(database='fireworks', collection_name='fireworks')
    out_dict, query_cursor = {}, []
    if args.workflows:
        workflows_db = D3Database(database='fireworks', collection_name='workflows')
        fw_ids = []
        if args.state:
            query_cursor = workflows_db.coll.find({"state": args.state}).distinct("fw_states")
            [fw_ids.extend(i.keys()) for i in query_cursor]
        fw_ids = list(set([int(i) for i in fw_ids]))
        query_cursor = fireworks_db.coll.find({"fw_id": {"$in": fw_ids}}).distinct("spec.identifier")
    elif args.special:
        query_cursor = frontend_db.coll.find({
            "mol_characterization.oxidation_potential": {'$exists': True},
            "species_characterization.cation1.spectra": {'$exists': True},
            "mol_characterization.reduction_potential": {'$exists': False},
            "species_characterization.anion1.geometry": {'$exists': False},
            # "mol_characterization.oxidation_potential": {'$exists': False},
            # "species_characterization.groundState.geometry": {'$exists': True},
            # "species_characterization.cation1.geometry": {'$exists': True},
            # '$or': [
            #     {"species_characterization.groundState.solvation_energy": {'$exists': False}},
            #     {"species_characterization.cation1.solvation_energy": {'$exists': False}},
            # ]
        }).distinct("_id")
    elif args.mol_characterization:
        search = "mol_characterization." + args.mol_characterization
        query_cursor = frontend_db.coll.find({search: {'$exists': False}}).distinct("_id")
    elif args.fw_special:
        query_cursor = fireworks_db.coll.find({}).distinct("spec.identifier")
        print(query_cursor)
    elif args.state:
        query_cursor = fireworks_db.coll.find({"state": args.state}).distinct("spec.identifier")
    elif args.unfinsihed:
        query_cursor = frontend_db.coll.find({'$or': [
            {"mol_characterization.oxidation_potential": {'$exists': False}},
            {"species_characterization.cation1.spectra": {'$exists': False}},
            {"species_characterization.cation2.spectra": {'$exists': False}},
            {"mol_characterization.reduction_potential": {'$exists': False}},
            {"species_characterization.anion1.spectra": {'$exists': False}},
            {"species_characterization.anion2.spectra": {'$exists': False}},
        ]
        }).distinct("_id")

    id_list = set([i for i in query_cursor])
    id_list = list(set(id_list))

    if args.not_in:
        with open(args.not_in) as fn:
            not_dict = json.load(fn)
            not_list = not_dict.keys()
        id_list = [i for i in id_list if i not in not_list]

    print(','.join(id_list))
    print("{} unfinished ids collected.".format(len(id_list)))
    title = "unfinished" if args.unfinsihed else (args.state or "").lower() or "collected"
    with open("reports/db_{}_ids_{}.json".format(title, time.strftime("%Y%m%d-%H%M%S")), 'w') as fn:
        json.dump({k: k for k in id_list}, fn, indent=2)
