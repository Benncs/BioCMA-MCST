from .handle_module import *
import os
import time
import numpy as np
from typing import Optional
import sys

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
    # from mpi4py import MPI

    # Set environment variables for OpenMP
    os.environ["OMP_PROC_BIND"] = "spread"
    os.environ["OMP_PLACES"] = "threads"
    # Initialize the BioMC handle
    handle = handle_module.init_handle(sys.argv)
    n_rank = handle_module.n_rank(handle)
    i_rank = handle_module.i_rank(handle)
    return handle, i_rank, n_rank


def init_simulation(
    outfolder: str,
    simulation_name: str,
    cma_path: str,
    params,
    is_recursive: bool = False,
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
    handle_module.register_cma_path(handle, cma_path, is_recursive)

    return handle


def set_initial_concentrations(
    handle, liquid: np.ndarray, gas: Optional[np.ndarray] = None
):
    if gas is not None:
        if liquid.shape != gas.shape:
            raise RuntimeError("Concentrations should be the same")

    handle_module.set_initialiser_from_data(handle, liquid.shape[0], liquid, gas)


def fast_run(
    outfolder: str,
    name: str,
    params: dict,
    cma_path: str,
    *,
    model_name: str,
    n_compartment: int,
    s_feed: float,
    liquid_flow_rate: float = 0.0,
    is_recursive: bool = False,
    is_serde: bool = False,
    f_init=None,
    serde_path=None,
):
    """Run the simulation with optional recursion or serde handling."""

    # Prepare parameters
    params = handle_module.make_params(**params)
    handle = handle_module.init_simulation(
        outfolder, name, cma_path, params, is_recursive
    )

    # Liquid feed configuration
    if liquid_flow_rate != 0:
        handle_module.set_liquid_feed_constant(handle, liquid_flow_rate, s_feed, 0, 0)

    # Model setup
    if not is_serde:
        if f_init is None:
            raise ValueError("f_init must be provided when is_serde is False")
        handle_module.set_initial_concentrations(handle, *f_init(n_compartment))
        handle_module.register_model_name(handle, model_name)
    else:
        if not serde_path:
            raise ValueError("serde_path must be provided when is_serde is True")
        handle_module.register_serde(handle, serde_path)

    # Apply and execute
    ok, msg = handle_module.apply(handle, is_serde)
    if not ok:
        print(msg)
        return -1

    handle_module.exec(handle)
    return 0


__all__.extend(
    ["pyinit_handle", "init_simulation", "set_initial_concentrations", "fast_run"]
)
