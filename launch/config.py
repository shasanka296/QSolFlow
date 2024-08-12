import os

config_solflow_file="CONFIG"
path_to_conda= input("Path to the miniconda directory (input as: path/miniconda_dir_name): ")
QSOl_flow_dir= input("Path to the QSolFlow dir (input as: path/to/solflow_dir): ")
is_this_HPC= input("Is this an HPC that has a singularity module? [y/n]")
if is_this_HPC.strip() =="y":
    loading_command= input("What is the loading command for singularity (i.e. module load ccs/singularity): ")
if is_this_HPC.strip() =="n":
    print("Please be sure that singularity is installed")
out_put_dir=  input("Path to where you want the outputs to be written: ")
Boss_dir = input("Path to boss dir: ")
singularity_path= input("Where is the f.sif file: ")
MDP_or_not=input("Do you have your MDP files? [y/n]")
MDP_path=None
Automatic=input("Do you want to automatically run the system after it is run[y/n]")
if MDP_or_not.strip() =="y":
    MDP_path=input("What is the path: ")
if MDP_or_not.strip()=="n":
    print("Go to {add_link} to see how to make one.")

monog_host= input("Paste the server host replacing in your password and adding the name of your DB: ")

mdp_path=os.path.join(MDP_path,"MDP")
source_conda_path=os.path.join(path_to_conda,"bin","activate")
envlines = [
    f'meta_dir="{out_put_dir}"\n',
    f'BOSSDIR="{Boss_dir}"\n',
    f'CONDAPATH="{source_conda_path}"\n',
    f'SINGPATH="{os.path.join(singularity_path,"f.sif")}"\n',
    f'MDP_Location= "{None if MDP_or_not=="n" else mdp_path}"\n'
]
to_write = [f'source {os.path.join(path_to_conda, "bin", "activate")}\n',
            'conda activate ASMD\n',
            'module load ccs/singularity\n',
            'module load gnu/5.4.0\n',
            'module load openmpi/1.10.7\n',
            'module load ccs/gromacs/skylake-gpu/2019\n',
            'source /opt/ohpc/pub/libs/gnu/openmpi/ccs/gromacs/2019/skylake-gpu/bin/GMXRC\n',
            f'export PYTHONPATH={os.path.join(QSOl_flow_dir,"QSolFlow")}:$PYTHONPATH']

Fw_config_lines = [
    "ECHO_TEST: FW_config.yaml is loaded for QSolFlow\n",
    "\n",
    "#ADD_USER_PACKAGES:\n",
    "#  -workflows.python\n",
    "\n",
    f"LAUNCHPAD_LOC: {os.path.join(QSOl_flow_dir,'QSolFlow','launch','md_launchpad.yaml')}\n",
    "\n",
    f"FWORKER_LOC: {os.path.join(QSOl_flow_dir,'QSolFlow','launch','my_fireworker.yaml')}\n",
    "\n",
    "\n",
    "# set to FIFO if you want older FireWorks to be run first, FILO if you want recent FireWorks run first.\n",
    "# Note that higher priority FireWorks are always run first.\n",
    "SORT_FWS: 'FIFO'\n",
    "\n",
    "# whether to print the FW.json file in your run directory\n",
    "PRINT_FW_JSON: True\n",
    "\n",
    "# whether to print the FW.yaml file in your run directory\n",
    "PRINT_FW_YAML: False\n",
    "\n",
    "# the name to give the script for submitting PBS/SLURM/etc. queue jobs\n",
    "SUBMIT_SCRIPT_NAME: FW_submit.sh\n",
    "\n",
    "# format for loggers (this String will be passed to logging.Formatter())\n",
    'FW_LOGGING_FORMAT: "%(asctime)s %(levelname)s %(message)s"\n',
    "\n",
    "# set True if you want the Queue Launcher to always create a new block directory every time it is called,\n",
    "# False if you want to re-use previous blocks\n",
    "ALWAYS_CREATE_NEW_BLOCK: False\n",
    "\n",
    "# where to store templates if you are using the TemplateWriterTask.\n",
    "# TEMPLATE_DIR:\n",
    "\n",
    "# tries to delete empty launch directories created when setting the _launch_dir in the spec of your Firework.\n",
    "REMOVE_USELESS_DIRS: True\n",
    "\n",
    "# if True, when rerunning a FIZZLED Firework, the serialized exception details are added to the spec.\n",
    "EXCEPT_DETAILS_ON_RERUN: True\n",
    "\n",
    "# the default host on which to run the web server\n",
    "WEBSERVER_HOST: 127.0.0.1\n",
    "\n",
    "# the default port on which to run the web server\n",
    "WEBSERVER_PORT: 5000\n",
    "#\n",
    "# the max length of the job name to send to the queuing system (some queuing systems limit the size of job names)\n",
    "QUEUE_JOBNAME_MAXLEN: 20\n"
]

launchpad = [
    f"host: {monog_host}\n",
    "logdir: null\n",
    "uri_mode: true\n"
]

my_firework = [
    "name: MD\n",
    "category: [\"gromacs\"]\n",
    "query: '{}'\n",
    "env:\n",
    f"    path: {os.path.join(QSOl_flow_dir,'launch')}\n",
    f"    runfile_log: {os.path.join(out_put_dir,'RunFile.log')}\n"
]




with open(f"{os.path.join(QSOl_flow_dir,'QSolFlow','solflow','workflows','envwf.py')}","w") as f:
    f.writelines(envlines)

with open(f"{os.path.join(QSOl_flow_dir,'QSolFlow','launch','FW_config.yaml')}","w") as f:
    f.writelines(Fw_config_lines)
with open(f"{os.path.join(QSOl_flow_dir,'QSolFlow','launch','md_launchpad.yaml')}","w") as f:
    f.writelines(launchpad)
with open(f"{os.path.join(QSOl_flow_dir,'QSolFlow','launch','my_fireworker.yaml')}","w") as f:
    f.writelines(my_firework)

with open(config_solflow_file,'w') as file:
    file.writelines(to_write)






desktop_file_path = os.path.join(QSOl_flow_dir, 'QSF.desktop')

source_file = f"{os.path.join(QSOl_flow_dir,'QSolFlow','launch','CONFIG')}"
python_script = f"{os.path.join(QSOl_flow_dir,'QSolFlow','solflow','Fast','fw_gui.py')}"
icon_path = f"{os.path.join(QSOl_flow_dir,'QSolFlow','launch','conda_env','QSF5.png')}"

command1 = f'bash -c "source {source_file} && python {python_script}"'
command2=  f'bash -c "source {source_file} && python {python_script} && cd QSolFlow/launch && rlaunch rapidfire"'

desktop_file_content = f"""[Desktop Entry]
Name=QSF
Comment=For help email sla296@uky.edu
Exec={command1 if Automatic.strip()=="n"else command2 }
Icon={icon_path}
Terminal=true
Type=Application
Categories=Utility
"""

with open(desktop_file_path, 'w') as file:
    file.write(desktop_file_content)

os.chmod(desktop_file_path, 0o755)

print(f'Desktop launcher created at {QSOl_flow_dir}')






