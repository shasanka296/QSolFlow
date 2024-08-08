import subprocess


class FIdel:
    def __init__(self, name, email,dir):
        self.dir = dir
        self.email=email
        self.name = name
        sumbit_name = self.name + 'Sumbit.sh'
        idel = subprocess.run(['sinfo'], shell=True, capture_output=True)
        nodes = idel.stdout.decode().splitlines()
        self.cacL_nodes = []

        file = [
            "#!/bin/bash",
            "",
            "",
            "",
            "#SBATCH -t 5-00:00:00                                   #Time for the job to run",
            f"#SBATCH --job-name={self.name}          #Name of the job",
            "",
            "#SBATCH -N 1                                    #Number of nodes required",
            "#SBATCH -n 40                           #Number of cores needed for the job",
            "#SBATCH --partition=CAC48M192_L         #Name of the queue",
            "",
            "#SBATCH --mail-type ALL                         #Send email on start/end",
            f"#SBATCH --mail-user {self.email}              #Where to send email",
            "",
            "#SBATCH --account=col_cmri235_uksr              #Name of account to run under",
            "",
            "#SBATCH --error=SLURM_JOB_%j.err                #Name of error file",
            "#SBATCH --output=SLURM_JOB_%j.out               #Name of output file",
            "",
            "#Module needed for this Gaussian job",
            f"#SBATCH --chdir={self.dir}/{self.name}script",
            "module load ccs/gaussian/g16-A.03/g16-haswell",
            "",
            "echo \"Job $SLURM_JOB_ID running on SLURM NODELIST: $SLURM_NODELIST\"",
            "",
            "#Gaussian Program execution command",
            "sleep 5",
            "g16 ept_p1_lpg.gjf > ept_p1_lpg.log"
        ]

        for a in range(len(nodes)):
            split = nodes[a].split()
            if split[4].strip() == "idle" and split[0].strip() == "CAC48M192_L":
                self.cacL_nodes.append(split[0])
            elif split[4].strip() == "idle" and split[0].strip() == "CAL48M192_L":
                self.cacL_nodes.append(split[0])
        try:
            self.idel_node = self.cacL_nodes[0]
        except:
            print("No idel node, making sumbit file using CAC48M192_L")
            self.idel_node = "CAC48M192_L"
            with open(f'{self.name}script/{sumbit_name}', 'a') as a:
                file[9] = '#SBATCH --partition=' + self.idel_node + '			#Name of the queue'
                file[27] = f'g16 {self.dir}/{self.name}script/{self.name}-final.gjf > {self.dir}/{self.name}script/{self.name}-final.log'
                for iteams in file:
                    a.write(iteams + '\n')

            print("Submit file is made, now submititng")
        else:
            print("found Idel node")
            self.idel_node = self.cacL_nodes[0]
            with open(f'{self.name}script/{sumbit_name}', 'a') as a:
                file[9] = '#SBATCH --partition=' + self.idel_node + '			#Name of the queue'
                file[27] = f'g16 {self.dir}/{self.name}script/{self.name}-final.gjf > {self.dir}/{self.name}script/{self.name}-final.log'
                for iteams in file:
                    a.write(iteams + '\n')

            print("Submit file is made, now submititng")
