from .handle_module import *
import os
import time
import numpy as np
from typing import Optional

__all__ = []
__doc__ = handle_module.__doc__
if hasattr(handle_module, "__all__"):
    __all__ = handle_module.__all__


def pyinit_handle(sim_id: int):
    """
    Initialize the BioMC handle.

    Returns:
        handle: The initialized BioMC handle.
    """
    from mpi4py import MPI

    # Set environment variables for OpenMP
    os.environ["OMP_PROC_BIND"] = "spread"
    os.environ["OMP_PLACES"] = "threads"

    # Get the MPI communicator
    comm = MPI.COMM_WORLD

    # Get the number of ranks and the rank ID
    n_rank = comm.Get_size()
    i_rank = comm.Get_rank()

    # Initialize the BioMC handle
    return handle_module.init_handle(n_rank, i_rank, sim_id, 1), i_rank, n_rank


def init_simulation(
    outfolder: str,
    simulation_name: str,
    cma_path: str,
    params,
    _id: Optional[int] = None,
):
    sim_id = 0
    if _id is None:
        sim_id = int(time.time())
    else:
        sim_id = _id

    full_out_dir = f"{outfolder}/{simulation_name}"

    os.makedirs(full_out_dir, exist_ok=True)

    # Initialize the BioMC handle
    handle, i_rank, n_rank = pyinit_handle(sim_id)

    # Register the simulation parameters
    handle_module.register_parameters(handle, params)

    # Register the result path
    handle_module.register_result_path(handle, f"{full_out_dir}/{simulation_name}")

    # Register the CMA path
    handle_module.register_cma_path(handle, cma_path)

    return handle


def set_initial_concentrations(
    handle, liquid: np.ndarray, gas: Optional[np.ndarray] = None
):
    if gas is not None:
        if liquid.shape != gas.shape:
            raise RuntimeError("Concentrations should be the same")

    handle_module.set_initialiser_from_data(handle, liquid.shape[1], liquid, gas)


__all__.extend(["pyinit_handle", "init_simulation","set_initial_concentrations"])
