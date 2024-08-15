import os
import subprocess
import time
class sumbit:
    def __init__(self, mole, userid, multicharge, dir):
        self.dir = dir
        self.mol = mole
        self.user = userid
        self.multi = multicharge


        while os.path.isfile(f'{self.dir}/{self.mol}script/{self.mol}.gjf') == False:
            time.sleep(1)
        with open(f'{self.dir}/{self.mol}script/{self.mol}.gjf', 'r+') as g, open(f'{self.dir}/{self.mol}script/{self.mol}-final.gjf', 'a') as p:
            a = g.readlines()
            b = 6
            l=[]       
            for lines in a:
                l.append(lines)

            D = [
                '%chk=EPT_opt_freq.chk\n',
                '%nprocshared=32\n',
                '%mem=4000MB\n',
                '# opt freq b3lyp/6-31g(d,p) scrf=(cpcm,solvent=acetonitrile) scf=qc\n',
                'scf=qc pop=(mk,full)\n'
                '\n',
                f'{self.mol}\n',
                '\n'


            ]
            try:
                c = [f'{int(self.multi)}', f' {int(self.multi) + 1}\n']
            except:
                c = ["0", " 1\n"]



            for lines in D:
                p.write(lines)
            for c_lines in (c):
                p.write(c_lines)
            for i in range(b, len(l)):
                p.write(l[i])
        subprocess.run(['sbatch', f'{self.dir}/{self.mol}script/{self.mol}Sumbit.sh'], check=True)
        a = subprocess.Popen(['squeue', '-u', f'{self.user}'], stdout=subprocess.PIPE)
        while True:
            output = a.communicate()[0].decode()
            if output.startswith("Priority") or output.startswith("Resources"):
                print(output)
            else:
                break
    
        print("dft is successfully submitted")
        while os.path.isfile(f'{self.dir}/{self.mol}script/{self.mol}-final.log') == False:
            time.sleep(1)
       	print("log is made, waiting for g16")

