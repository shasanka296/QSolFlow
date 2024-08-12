<img src="https://github.com/user-attachments/assets/69eac423-3c0e-4209-898f-1825ab612ac4" alt="QSF5" width="300"/>

## Setting Up and Running SolFlow

To get started with SolFlow, follow these steps:

1. **Download and Install Singularity**
   - If you're using a High-Performance Computing (HPC) environment that already has a Singularity module, you may skip this step.
   - For those needing to install Singularity, follow the instructions [here](https://docs.sylabs.io/guides/3.0/user-guide/installation.html).

2. **Set Up and Configure MongoDB**
   - Detailed instructions, including a video tutorial, can be found [here](#) (need to add link).

3. **Install Miniconda**
   - Download and install Miniconda by following the instructions provided [here](https://docs.anaconda.com/miniconda/miniconda-install/).

4. **Set Up the Conda Environments**
   - If you are using the classic solver, switch to the new libmamba solver using the following command:
     ```bash
     conda update -n base conda --solver=classic
     conda install -n base conda-libmamba-solver -y
     conda config --set solver libmamba
     ```
   - Use the following commands to create the necessary Conda environments:
     ```bash
     git clone https://github.com/shasanka296/QSolFlow.git
     conda env create -f QSolFlow/launch/conda_env/main_env.yml
     conda env create -f QSolFlow/launch/conda_env/DFT.yml
     conda env create -f QSolFlow/launch/conda_env/LIG.yml
     ```

   - Inside the QSolFLow directory download ligpargen:
     ```bash
     conda activate ligpg
     git clone https://github.com/Isra3l/ligpargen.git
     pip install -e ligpargen
     ```

6. **Run the Configuration Script**
   - Execute the configuration script located in the `launch` directory.


Follow these steps carefully to ensure SolFlow is set up and running correctly.

For any questions email sla296@uky.edu
