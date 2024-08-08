## Setting Up and Running SolFlow

To get started with SolFlow, follow these steps:

1. **Download and Install Singularity**
   - If you're using a High-Performance Computing (HPC) environment that already has a Singularity module, you may skip this step.
   - For those needing to install Singularity, follow the instructions [here](https://docs.sylabs.io/guides/3.0/user-guide/installation.html).

2. **Set Up and Configure MongoDB**
   - Detailed instructions, including a video tutorial, can be found [here](#) (add your video link).

3. **Install Miniconda**
   - Download and install Miniconda by following the instructions provided [here](https://docs.anaconda.com/miniconda/miniconda-install/).

4. **Set Up the Conda Environment**
   - Navigate to the `launch/conda_env` directory and set up the environment using the provided files.

5. **Run the Configuration Script**
   - Execute the configuration script located in the `launch` directory.

6. **Run SolFlow**
   - Ensure the `SolFlow.sh` script has the correct permissions:
     ```bash
     chmod 777 SolFlow.sh
     ```
   - Run the script:
     ```bash
     ./SolFlow.sh
     ```

Follow these steps carefully to ensure SolFlow is set up and running correctly.
For any questions email sla296@uky.edu
