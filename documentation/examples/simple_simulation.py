"""
Minimal working program to use BioMC API V0.5.

This program demonstrates a basic example of how to perform a simulation using BioMC from a Python program.
It focuses on:
- Creating a handle for shared execution (No MPI)
- Performing basic simulation settings and running it
- Loading and running a basic simulation

Note that the linked library is the shared one, this will not work linking the distributed one.

Usage:
    Run this script to perform a simulation using BioMC API V0.5.
"""

# Import the handle module
import handle_module
import numpy as np

# define the constants
OUTFOLDER = "./out/"
SIMULATION_NAME = "my_simulation_name"
CMA_PATH = "/path/to/the/cma/"  # don´t forget last /


def run(params):
    handle = handle_module.init_simulation(OUTFOLDER, SIMULATION_NAME, CMA_PATH, params)
    # cma with 500 compartment, simulation with 4 species liquid only
    liquid_concentration_0 = np.zeros((500, 4))
    handle_module.set_initial_concentrations(handle, liquid_concentration_0)
    handle_module.register_model_name(handle, "model_name")
    # Apply the simulation settings
    rc = handle_module.apply(handle, False)

    # Check if the simulation settings were applied successfully
    if not rc[0]:
        print(rc[1])
        return -1

    # Execute the simulation
    rc = handle_module.exec(handle)


if __name__ == "__main__":
    # Create the simulation parameters
    params = handle_module.make_params(
        biomass_initial_concentration=0,
        final_time=250,
        delta_time=0.5,
        number_particle=0,
        number_exported_result=5,
    )

    run(params)
