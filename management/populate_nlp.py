"""
Parse and push data from NLP parsing and submit calculations for extracted molecules
"""
import argparse
from monty.serialization import loadfn, dumpfn
from d3tales_api.D3database.d3database import *
from d3tales_api.Processors.back2front import *
from d3tales_api.Processors.d3tales_parser import *
from d3tales_fw.management.populate_lpad import populate_d3tales_lpad

parser = argparse.ArgumentParser(description='Parse and push data from NLP parsing and submit calculations for extracted molecules ')
parser.add_argument('filename', metavar='filename', type=str, help='filepath for a JSON nlp data file', default=False)
parser.add_argument('-m', '--nlp_model', type=str, help='NLP model name', default="test_model")
parser.add_argument('-p', '--priority', type=int, help='jobs priority', default=5)
parser.add_argument('-nc', '--no_calculations', action='store_true', help='no calculations if used', default=5)
args = parser.parse_args()

all_nlp_data = loadfn(args.filename)
all_ids = {}
for doi, html_str in all_nlp_data.items():
    # Process html string and produce backend data
    instance = ProcessNlp.from_html(html_str, args.nlp_model, doi=doi).data_dict
    # Inserts backend data
    back_id = BackDB(collection_name='nlp', instance=instance, last_updated=True).id
    # Process backend data to frontend data and insert frontend data
    mol_ids = DOI2Front(doi=back_id, insert=True).mol_ids

    # Initiate high-throughput DFT calculations
    if args.no_calculations:
        fw_ids = ["" for _ in mol_ids]
    else:
        fw_ids = populate_d3tales_lpad(ids_list=mol_ids, priority=args.priority)

    all_ids.update(zip(mol_ids, fw_ids))

dumpfn(all_ids, "nlp_data/nlp_{}_molIDs.json".format(args.filename.split("/")[-1]), indent=2)
