import tkinter as tk
from tkinter import ttk
import scriptMaker as m
import lig as l
import dft as dft
import waiting as wait
import chargeTrasnfer as transfer
import packmol as pack
import make_gro as gro
import subprocess
import ASMD_1 as ASMD


class GUI:
    def __init__(self):
        self.window = tk.Tk()

        self.window.title('FastAtom')

        self.window.geometry('850x650')

        self.button = tk.Button(self.window,
                                text="Submit",
                                width=5,
                                height=1,
                                bg="gray",
                                fg="black", command=self.run)

        a= subprocess.run("echo $PWD", shell=True, capture_output=True)
        self.curentDirectory = a.stdout.decode().strip("\n")
        self.SoluteSmiles = tk.Entry(self.window, fg="black", bg="white", width=50)

        self.Density = tk.Entry(self.window,fg="black", bg="white", width=50 )
        self.xdim =tk.Entry(self.window,fg="black", bg="white", width=50 )
        self.ydim=tk.Entry(self.window,fg="black", bg="white", width=50 )
        self.zdim=tk.Entry(self.window,fg="black", bg="white", width=50 )

        self.xdimLabel= tk.Label(self.window, text="x dimensions: ".ljust(20))
        self.ydimLabel = tk.Label(self.window, text="y dimensions: ".ljust(20))
        self.zdimLabel = tk.Label(self.window, text="z dimenisons: ".ljust(20))

        self.email= tk.Entry(self.window, fg="black", bg="white", width=50)

        self.emailLable= tk.Label(self.window, text="Email: ".ljust(20))

        self.labelDensity= tk.Label(self.window, text ="Desity of solvent(M)".ljust(20))

        self.SoluteName = tk.Entry(self.window, fg="black", bg="white", width=50)


        self.SolventSmiles = tk.Entry(self.window, fg="black", bg="white", width=50)

        self.Solname= tk.Entry(self.window, fg= "black", bg="white", width=50)

        self.Concentrations = tk.Entry(self.window, fg="black", bg="white", width=50)

        self.labelSolute = tk.Label(self.window, text="Solute names (3 letters, if using multiple, seprate with commas): ".ljust(20))

        self.labelSolvent= tk.Label(self.window, text="Solvent name (3 letters): ".ljust(20))


        self.solventname2 = tk.Entry(self.window, fg="black", bg="white", width=50)

        self.solvent2Lable = tk.Label(self.window, text="Solvent2 name (3 letters): ")

        self.solventsmiles2 = tk.Entry(self.window, fg="black", bg="white", width=50)

        self.solvent2SmilesLable = tk.Label(self.window, text="Solvent2 Smiles: ")


        self.Density2= tk.Entry(self.window, fg="black", bg="white", width=50)
        self.Density2label=tk.Label(self.window, text ="Density of Solvent2(M): ".ljust(20))


        self.labelConcentration= tk.Label(self.window, text ="Concentration(M) of the Solutes, seprate with commas: ".ljust(20))

        self.lableSolventSmiles = tk.Label(self.window, text="Solvent SMILES code:".ljust(20))

        self.labelSoluteSmiles = tk.Label(self.window, text= "Solute SMILES code(separate with commas): ".ljust(20))

        self.user = tk.Entry(self.window, fg="black", bg="white", width=50)

        self.user_label = tk.Label(self.window, text="User ID")


        self.chargeV = tk.IntVar()

        self.ratio= tk.Entry(self.window, fg="black", bg="white", width=50)
        self.rationlabel= tk.Label(self.window, text="Ratio of the Solvent1/Solvent2")
        self.molarMass= tk.Entry(self.window, fg="black", bg="white", width=50)
        self.molarlabel= tk.Label(self.window, text="Molarmass of Solvents")
    

        self.chargeCheck = tk.Checkbutton(self.window, text="Charge on the Solutes", variable=self.chargeV)
        self.chargeMatrix = tk.Entry(self.window, fg="black", bg="white", width=50)

        self.chargeEntry = tk.Entry(self.window, fg="black", bg="white", width=50)

        self.dftV = tk.IntVar()

        self.DFTCheck = tk.Checkbutton(self.window, text="DFT", variable=self.dftV)

        self.progressBar= ttk.Progressbar(self.window, length= 400 , mode="determinate")
        self.progresslabel = tk.Label(self.window, text="Progress:")
        self.step = tk.Label(self.window, text="Not started",font=("Arial", 25))

        self.running = ttk.Progressbar(self.window, length= 400 , mode="indeterminate")


        self.labelSolvent.grid(row=0, column=0)
        self.Solname.grid(row=0, column=1)

        self.lableSolventSmiles.grid(row=1, column=0)
        self.SolventSmiles.grid(row=1, column=1)

        self.labelSolute.grid(row=2, column=0)
        self.SoluteName.grid(row=2, column=1)

        self.labelSoluteSmiles.grid(row=3, column=0)
        self.SoluteSmiles.grid(row=3, column=1)

        self.labelConcentration.grid(row=4, column=0)
        self.Concentrations.grid(row=4, column=1)

        self.user_label.grid(row=5, column=0)
        self.user.grid(row=5, column=1)

        self.labelDensity.grid(row=6, column=0)
        self.Density.grid(row=6, column=1)

        self.xdimLabel.grid(row=9, column=0)
        self.xdim.grid(row=9, column=1)

        self.ydimLabel.grid(row=10, column=0)
        self.ydim.grid(row=10, column=1)

        self.zdimLabel.grid(row=11, column=0)
        self.zdim.grid(row=11, column=1)

        self.chargeCheck.grid(row=12, column=0)
        self.chargeMatrix.grid(row=12, column=1)

        self.solvent2Lable.grid(row=13, column=0)
        self.solventname2.grid(row=13, column=1)

        self.solvent2SmilesLable.grid(row=14, column=0)
        self.solventsmiles2.grid(row=14, column=1)


        self.Density2.grid(row=15, column=1)
        self.Density2label.grid(row=15, column=0)

        self.rationlabel.grid(row=16, column=0)
        self.ratio.grid(row=16, column=1)

        self.emailLable.grid(row=20, column=0)
        self.email.grid(row=20, column=1)

        self.molarMass.grid(row=21, column=1)
        self.molarlabel.grid(row=21, column=0)
    

        self.DFTCheck.grid(row=22, column=0)

        self.button.grid(row=23, column=0)

        self.progressBar.grid(row=24, column=1)
        self.progresslabel .grid(row=24, column=0)


        self.step.grid(row=26, column=1)
        


        self.window.mainloop()

    def update_progress(self, txt, int):
        self.progressBar['value'] += int
        self.step.config(text=txt)
        self.window.update()
    def starting(self):
        self.running.start(100)
    def stop(self):
        self.running.stop()
    def run(self):
        self.starting()
        self.update_progress("Program started", 0)
        
        SoluteMatrix = self.SoluteName.get().split(",")
        SoluteSmilesMatrix = self.SoluteSmiles.get().split(",")
        charge=self.chargeMatrix.get().strip(",")

        for I in range(len(SoluteMatrix)):
            self.update_progress(f"Starting LigPG for solute {I+1}", 0)
            l.lig(SoluteSmilesMatrix[I], SoluteMatrix[I].strip()[:3] + f"_Solute{I+1}", charge[i],self.curentDirectory)

        self.update_progress(f"Starting LigPG for solvent1", 0)
        l.lig(self.SolventSmiles.get(), self.Solname.get()[:3] +"_Solvent",0, self.curentDirectory)

        if len(self.solventsmiles2.get()) != 0:
            self.update_progress(f"Starting LigPG for solvent2", 0)
            l.lig(self.solventsmiles2.get(), self.solventname2.get()[:3] + "_Solvent2", 0,self.curentDirectory)

        if self.dftV.get()==1:
            print("dft was selcted")
            self.update_progress("Starting DFT",25)
            self.script(SoluteMatrix)

        else:
            print("dft not selcted")
            self.update_progress("DFT not selted, starting packmol", 25)
            self.pack(SoluteMatrix)
    def script(self, solutename):
        SoluteM = solutename
        charge = self.chargeMatrix.get().strip(",")
        self.update_progress("Starting solvent 1 DFT", 10)

        m.FIdel(self.Solname.get()[:3] +"_Solvent", self.email.get(), self.curentDirectory)
        self.update_progress("Script made for solvent 1",0)
        dft.sumbit(self.Solname.get()[:3] +"_Solvent", self.user.get(), "0", self.curentDirectory)
        self.update_progress("Dft sumbited for Solvent 1, wating on g16",0)
        wait.waiting(self.Solname.get()[:3] +"_Solvent", self.curentDirectory)
        self.update_progress ("Transfering charges for Solvetn 1",0)
        transfer.trans(self.Solname.get()[:3] +"_Solvent", self.curentDirectory)
        self.update_progress ("charge transfer done for Solvent 1, starting the DFT for the solutes",15)



        for i in range(len(SoluteM)):
            m.FIdel(SoluteM[i].strip()[:3] + f"_Solute{i+1}", self.email.get(), self.curentDirectory)
            self.update_progress(f"Script made for solute{i+1}",0)
            dft.sumbit(SoluteM[i].strip()[:3] + f"_Solute{i+1}", self.user.get(), charge.strip(),
                       self.curentDirectory)
            self.update_progress(f"Dft sumbited for Solute {i+1}, wating on g16",0)
            wait.waiting(SoluteM[i].strip()[:3] + f"_Solute{i+1}", self.curentDirectory)
            self.update_progress(f"Dft sumbited for Solute{i+1}, wating on g16",0)
            transfer.trans(SoluteM[i].strip()[:3] + f"_Solute{i+1}", self.curentDirectory)
            self.update_progress(f"Dft sumbited for Solute {i+1}, wating on g16",0)


        if len(self.solventsmiles2.get()) != 0:
            m.FIdel(self.solventname2.get()[:3] + "_Solvent2", self.email.get(), self.curentDirectory)
            self.update_progress("Script made for solvent 2", 0)
            dft.sumbit(self.solventname2.get()[:3] + "_Solvent2", self.user.get(), "0",
                       self.curentDirectory)
            self.update_progress("Dft sumbited for Solvent 2, wating on g16",0)
            wait.waiting(self.solventname2.get()[:3] + "_Solvent2", self.curentDirectory)
            self.update_progress("Transfering charges for Solvetn 2",0)
            transfer.trans(self.solventname2.get()[:3] + "_Solvent2", self.curentDirectory)
            self.update_progress("charge transfer done for Solvent 2, now using packmol",25)

        self.pack(SoluteM)
    def pack(self, solutename):
        solutes = solutename
        SoluteConetrationMatrix = self.Concentrations.get().split(",")
        pack.Solvate(self.Solname.get()[:3], solutes, SoluteConetrationMatrix, self.Density.get(), self.solventname2.get()[:3], self.Density2.get(), self.xdim.get(), self.ydim.get(), self.zdim.get(), self.curentDirectory, self.ratio.get())
        self.update_progress("Packmol done, creating gro file and cleaing up",10)
        gro.gro(self.Solname.get()[:3], solutes,self.solventname2.get()[:3], self.curentDirectory, self.xdim.get(),self.ydim.get(),self.zdim.get())
        self.update_progress("starting simulation",0)
        runer=ASMD.ASMD()
        a=runer.EnergyMin()
        self.update_progress(f"Energy Minimization done, exited wiht:{a}",0)
        b=runer.NVT()
        self.update_progress(f"NVT step done exited with:{b}",0)
        c=runer.NPT()
        self.update_progress(f"NPT done exited with:{c}",0)
        Denisty=runer.calculate_density()
        self.update_progress(f"Densities: {Denisty}",0)
        d=runer.check_density_accuracy( float(self.Density.get())* float(self.molarMass.get()),Denisty)
        #d=runer.check_density_accuracy( 0.33,Denisty)
        self.update_progress(f"checking density, exited with:{d}",0)
        e=runer.correct()
        self.update_progress(f"PBC correction complete: e",0)
        f=runer.index_file()
        self.update_progress(f"Index file made exited with:{f}",0)
        g=runer.msd_mkaer()
        self.update_progress(f"MSD file made exited wiht:{g}",0)
        residues=runer.extract_residues_from_itp()
        self.update_progress(f"Residue names: {residues}",0)
        cord=runer.rdf(residues,5)
        self.update_progress(f"cordiantion numer{cord}",0)
        h=runer.cordination_number(cord)
        self.update_progress(f"Program finished exited wiht:{h} ", 100)
        self.stop()


GUI()



