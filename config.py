import os
config_solflow_file="CONFIG"
path_to_conda= input("Path to the miniconda directory: ")
SOl_flow_dir= input("Path to the SolFlow dir: ")
is_this_HPC= input("Is this an HPC that has a singularity module? [y/n]")
if is_this_HPC =="y":
    loading_command= input("what is the loading command for singularity (i.e. module load ccs/singularity): ")
if is_this_HPC =="n":
    print("Please be sure that singularity is installed, SolFlow configuration done")

to_write=[f'source {os.path.join(path_to_conda,"miniconda3","bin","activate")}\n',
'conda activate ASMD\n',
'module load ccs/singularity\n',
'module load gnu/5.4.0\n',
'module load openmpi/1.10.7\n',
'module load ccs/gromacs/skylake-gpu/2019\n',
'source /opt/ohpc/pub/libs/gnu/openmpi/ccs/gromacs/2019/skylake-gpu/bin/GMXRC\n',
f'export PYTHONPATH={SOl_flow_dir}:$PYTHONPATH']

with open(config_solflow_file,'a') as file:
    file.writelines(to_write)

