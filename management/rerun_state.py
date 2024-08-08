import os
from tqdm import tqdm
import multiprocessing
from pathlib import Path
from fireworks import LaunchPad
BASE_DIR = Path(__file__).resolve().parent.parent

STATE = "FIZZLED"
FW_IDS = ""
TIMEOUT = 10  # timout for a rerun request in seconds

lpad_file = os.path.join(BASE_DIR.parent, 'config', 'my_launchpad.yaml')
lpad = LaunchPad().from_file(lpad_file)

fw_ids = [int(i) for i in FW_IDS.split(",")] if FW_IDS else lpad.get_fw_ids({"state": STATE})
for fw_id in tqdm(fw_ids):
    p = multiprocessing.Process(target=lpad.rerun_fw, name=f"rerun_fw_{fw_id}", args=(fw_id,))
    p.start()
    p.join(TIMEOUT)

    # If thread is active
    if p.is_alive():
        print(f"rerun_fw_{fw_id} is still running after {TIMEOUT} seconds... let's kill it...")
        # Terminate foo
        p.terminate()
        p.join()

