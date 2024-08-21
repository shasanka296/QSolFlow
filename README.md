<div style="text-align: center;">
   <img src="https://github.com/user-attachments/assets/693cd197-f5df-48bb-9741-17f37b7aeacc" alt="QSF_web" width="400">
</div>

# QSolFlow: Charged (Q) Solution (Sol) Dynamics Workflow (Flow) Automator

## Setting Up and Running QSolFlow

To get started with QSolFlow, follow these steps:

***A Video guide to assist with the download process and setup can be found [here](https://www.dropbox.com/scl/fi/5anxgmcrh07idj4h1w9yi/Download-QSolFLow.mp4?rlkey=tah2rtelcvm5bqowcbgb6sa87&st=kyiy4w4n&dl=0)***

1. **Download and Install Singularity**
   - If you're using a high-performance computing (HPC) environment that already has a Singularity module, you may skip this step.
   - For those needing to install Singularity, follow the instructions [here](https://docs.sylabs.io/guides/3.0/user-guide/installation.html).

2. **Set Up and Configure MongoDB, plus video instructions**
   - Detailed instructions, including a video tutorial, to set up MongoDB can be found [here](https://www.dropbox.com/scl/fi/pw6toktp4jyqye7qqr684/Database.mp4?rlkey=wh1ecsenzzqpruppruwksrren&st=1mev4rpr&dl=0).
     
3. **Install Miniconda**
   - Download and install Miniconda by following the instructions provided [here](https://docs.anaconda.com/miniconda/miniconda-install/).

4. **Set Up the Conda Environments**
   - If you are using the classic solver, switch to the new libmamba solver using the following command:
     ```bash
     conda update -n base conda --solver=classic
     y
     conda install -n base conda-libmamba-solver
     y
     conda config --set solver libmamba
     
     ```
   - Use the following commands to create the necessary conda environments:
     ```bash
     conda install git
     git clone https://github.com/shasanka296/QSolFlow.git
     conda env create -f QSolFlow/launch/conda_env/main_env.yml
     conda env create -f QSolFlow/launch/conda_env/DFT.yml
     conda env create -f QSolFlow/launch/conda_env/LIG.yml
     
     ```

   - Inside the QSolFLow directory download ligpargen:
     ```bash
     cd QSolFlow
     conda activate ASMD
     git clone https://github.com/Isra3l/ligpargen.git
     conda activate ligpg
     pip install -e ligpargen
     
     ```

6. **Run the Configuration Script**
   - When setting the output directory make sure that the total path is less than 80 characters

   - Execute the configuration script located in the `launch` directory.
      ```bash
     conda activate base
     cd launch
     python config.py
     
     ```
8. **Download the singularity container, MDP files and boss through here along with some example files***
   - [DropBox link](https://www.dropbox.com/scl/fo/bcmwn6ufjk6k5qrt58s36/AC3_o6bwXv0xvVED3PitmX0?rlkey=mt82tc7ampn2ts1tui6gsx5ti&st=flshzr1w&dl=0)


Follow these steps carefully to ensure SolFlow is set up and running correctly.

For any questions email sla296@uky.edu or antonsperea@uky.edu

