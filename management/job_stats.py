import argparse
from tqdm import tqdm
from d3tales_api.D3database.d3database import D3Database


def get_runtimes(search_term, mcc_v_stmp=True):
    fireworks_db = D3Database(database='fireworks', collection_name='fireworks')
    launch_db = D3Database(database='fireworks', collection_name='launches')

    fireworks_cursor = fireworks_db.coll.find({'$and': [
        {"state": "COMPLETED"},
        {"name": {'$regex': search_term}}
    ]})

    launch_ids = []
    for fw in fireworks_cursor:
        launch_ids.extend(fw.get("launches", []))

    stampede_times = []
    mcc_times = []
    for launch_id in tqdm(launch_ids):
        launch_cursor = launch_db.coll.find_one({"launch_id": launch_id})
        if launch_cursor.get("ip").startswith("206.76"):
            stampede_times.append(launch_cursor.get("runtime_secs"))
        elif launch_cursor.get("ip").startswith("10.30"):
            mcc_times.append(launch_cursor.get("runtime_secs"))

    if mcc_v_stmp:
        print(search_term.capitalize(), " Calculations: ")
        avg_mcc_time = round(sum(mcc_times) / len(mcc_times), 3)
        print("Total Average Runtime for {} jobs on MCC: \t\t{} min  ({} hr)".format(len(mcc_times), round(avg_mcc_time / 60, 3),
                                                                         round(avg_mcc_time / 3600, 3)))
        avg_stampede_time = round(sum(stampede_times) / len(stampede_times), 3)
        print("Total Average Runtime for {} jobs on Stampede: \t{} min  ({} hr)".format(len(stampede_times), round(avg_stampede_time / 60, 3),
                                                                            round(avg_stampede_time / 3600, 3)))
    total_avg_time = round(sum(mcc_times + stampede_times) / len(mcc_times + stampede_times), 3)
    print("+--------------------------------------------------------------------------+")
    print("|   Average {} Runtime for {} jobs: \t{} min  ({} hr) \t|".format(search_term.capitalize(), len(mcc_times+stampede_times), round(total_avg_time / 60, 3),
                                                                 round(total_avg_time / 3600, 3)))
    print("+--------------------------------------------------------------------------+")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Get stats for D3TaLES FireWorks workflow jobs.')
    parser.add_argument('search_term', metavar='search_term', type=str, help='type of terms jobs to get stats for')
    parser.add_argument('-t', '--totals', action='store_false', help='do not separate analysis by MCC vs Stampede2')
    args = parser.parse_args()

    if args.search_term == 'all':
        for term in ["opt", "energy", "tddft", "wtuning"]:
            get_runtimes(term, mcc_v_stmp=args.totals)
    else:
        get_runtimes(args.search_term)
