import matplotlib.pyplot as plt
from QSFlow.workflows.envwf import meta_dir
import os
def TitrationPlotter( Average_simulation_density_dic, key):
    """
     Plots the titration graph and saves it to the output dir
        Parameters:
            - Average_simulation_density_dic: Dictionary of average simulation density.
            - key: The key of the current system.
    """
    difrence=[float(f'{float(i):.4f}') for i in Average_simulation_density_dic.values() ]
    x_axis=[i for i in Average_simulation_density_dic]
    print(difrence)

    plt.plot(x_axis,difrence,'o')
    plt.xlabel("System")
    plt.ylabel("Simulation density")
    plt.title(f"Titration Result{key}")
    plt.savefig( os.path.join(meta_dir,f"Output{key}" f'TitrationGraph_{key}.pdf'))
    plt.clf()

