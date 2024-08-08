from tqdm import tqdm
from d3tales_api.D3database.d3database import D3Database

SPECIAL_KW = ""  # PUGREST.ServerBusy not successfully submitted ValidationError FileNotFound RESTAPI
RUN_MCC = False
RERUN = False or RUN_MCC or SPECIAL_KW
MULTIPLICITY = None
LAUNCH_CUTOFF = 0
PRINT_IDS = False
FW_STATE = "FIZZLED"

ERROR_KW1 = SPECIAL_KW or "not terminated"
ERROR_KW2 = SPECIAL_KW or "Structure not converged"
ERROR_KW3 = SPECIAL_KW or "Frequency calculation failed Normal termination"
ERROR_KW4 = SPECIAL_KW or "Calculation has a spin state greater than triplet"
ERROR_KW5 = SPECIAL_KW or "not possible for this molecule"
ERROR_KW6 = SPECIAL_KW or "molecule has no atoms"

if __name__ == "__main__":
    fizzled_ids, mol_ids = [], []
    fireworks_db = D3Database(database='fireworks', collection_name='fireworks')
    launch_db = D3Database(database='fireworks', collection_name='launches')

    if LAUNCH_CUTOFF:
        fizzled_fws = fireworks_db.coll.find({"state": "FIZZLED"})
        fizzled_fw_ids = [f.get("fw_id") for f in fizzled_fws]

        for fw_id in fizzled_fw_ids:
            if launch_db.coll.count_document({"fw_id": fw_id}) > LAUNCH_CUTOFF:
                fizzled_ids.append(fw_id)

    else:
        if FW_STATE == "FIZZLED":
            query_txt = {"state": FW_STATE,
                         '$or': [
                             # {}
                             {"action.stored_data._exception._stacktrace": {'$regex': ERROR_KW1}},
                             {"action.stored_data._exception._stacktrace": {'$regex': ERROR_KW2}},
                             {"action.stored_data._exception._stacktrace": {'$regex': ERROR_KW3}},
                             {"action.stored_data._exception._stacktrace": {'$regex': ERROR_KW4}},
                             {"action.stored_data._exception._stacktrace": {'$regex': ERROR_KW5}},
                             {"action.stored_data._exception._stacktrace": {'$regex': ERROR_KW6}},
                         ]
                         }
        else:
            query_txt = {"state": FW_STATE}
        launch_cursor = launch_db.coll.find(query_txt)

        uncoverged_fw_ids = set([l.get("fw_id") for l in launch_cursor])

        for fw_id in tqdm(uncoverged_fw_ids):
            if not fireworks_db.coll.count_documents({"fw_id": fw_id}):
                continue
            if fireworks_db.coll.find({"fw_id": fw_id})[0].get("state") == FW_STATE:
                fizzled_ids.append(fw_id)
                if RUN_MCC:
                    fireworks_db.coll.update_one({"fw_id": fw_id}, {"$set": {"spec._category": "gaussian_mcc"}})
                if RERUN:
                    fireworks_db.coll.update_one({"fw_id": fw_id},
                                                 {"$set": {"spec._tasks.0.check_if_already_run": False}})
                if MULTIPLICITY:
                    fireworks_db.coll.update_one({"fw_id": fw_id},
                                                 {"$set": {"spec._tasks.0.check_if_already_run": False}})
                if PRINT_IDS:
                    mol_ids.append(
                        fireworks_db.coll.find_one({"fw_id": fw_id}, {"spec.identifier": True}).get("spec", {}).get(
                            "identifier"))
    fizzled_ids = set(fizzled_ids)
    if RERUN:
        print('lpad rerun_fws -i ', ','.join(map(str, fizzled_ids)))
    else:
        print('lpad defuse_fws -i ', ','.join(map(str, fizzled_ids)))

    if PRINT_IDS:
        print(','.join(map(str, fizzled_ids)))

    print("Found {} unconverged firworks.".format(len(fizzled_ids)))
