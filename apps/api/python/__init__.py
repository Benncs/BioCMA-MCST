import os
import sys
from collections.abc import Callable
from typing import Optional

import numpy as np

from .handle_module import *

__all__ = []
__doc__ = handle_module.__doc__
if hasattr(handle_module, "__all__"):
    __all__ = handle_module.__all__

__all__.extend(
    [
        "pyinit_handle",
        "init_simulation",
        "set_initial_concentrations",
        "fast_run",
        "config_and_run",
        "check_version",
    ]
)


def check_version():
    if handle_module.get_version() != [1, 2, 1]:
        print(
            f"Warning: Script was written for v1.2.1, used {handle_module.get_version()}"
        )


def set_sim_env(num_threads: str = "1", udf_path: str = None, **kwargs):
    os.environ["KOKKOS_NUM_THREADS"] = num_threads
    if udf_path is not None:
        os.environ["BIOMC_LIB_UDF"] = udf_path
    for k, v in kwargs.items():
        os.environ[k] = v


def pyinit_handle(sim_id: int):
    """
    Initialize the BioMC handle.

    Returns:
        handle: The initialized BioMC handle.
    """

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
    sim_id: Optional[int] = None,
):
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

    handle_module.set_initialiser_from_data(handle, liquid.shape[0], liquid, gas)


def config_and_run(
    outfolder: str,
    name: str,
    params: dict,
    cma_path: str,
    model_name: str,
    is_serde: bool,
    callback_pre: Optional[Callable] = None,
    callback_post: Optional[Callable] = None,
):
    if callback_pre is not None:
        callback_pre(outfolder, name, params, model_name)

    # Prepare parameters
    params = handle_module.make_params(**params)

    handle = init_simulation(outfolder, name, cma_path, params)

    handle_module.register_model_name(handle, model_name)

    if callback_post is not None:
        callback_post(handle)

    # Apply and execute
    ok, msg = handle_module.apply(handle, is_serde)
    if not ok:
        print(msg)
        return -1

    handle_module.exec(handle)
    return handle_module.i_rank(handle)


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
    is_serde: bool = False,
    f_init=None,
    serde_path=None,
    uniform=None,
    extra=None,
):
    """Run the simulation with optional recursion or serde handling."""

    def cb_pre(_outfolder, _name, _params, _model_name):
        if uniform is not None:
            _params.uniform_particle_init = True

    def cb_post(handle):
        # Liquid feed configuration
        if liquid_flow_rate != 0:
            handle_module.set_liquid_feed_constant(
                handle, liquid_flow_rate, s_feed, 0, 0
            )
        # Model setup
        if not is_serde:
            if f_init is None:
                raise ValueError("f_init must be provided when is_serde is False")
            set_initial_concentrations(handle, *f_init(n_compartment))
        else:
            if not serde_path:
                raise ValueError("serde_path must be provided when is_serde is True")
            handle_module.register_serde(handle, serde_path)

        if extra is not None:
            extra(handle)

    return config_and_run(
        outfolder, name, params, cma_path, model_name, is_serde, cb_pre, cb_post
    )


# def fast_run(
#     outfolder: str,
#     name: str,
#     params: dict,
#     cma_path: str,
#     *,
#     model_name: str,
#     n_compartment: int,
#     s_feed: float,
#     liquid_flow_rate: float = 0.0,
#     is_serde: bool = False,
#     f_init=None,
#     serde_path=None,
#     uniform=None,
# ):
#     """Run the simulation with optional recursion or serde handling."""

#     # Prepare parameters
#     params = handle_module.make_params(**params)

#     if uniform is not None:
#         params.uniform_particle_init = True

#     handle = init_simulation(outfolder, name, cma_path, params)

#     # Liquid feed configuration
#     if liquid_flow_rate != 0:
#         handle_module.set_liquid_feed_constant(handle, liquid_flow_rate, s_feed, 0, 0)

#     handle_module.register_model_name(
#         handle, model_name
#     )  # Needed to set it even if serde with UDF

#     # Model setup
#     if not is_serde:
#         if f_init is None:
#             raise ValueError("f_init must be provided when is_serde is False")
#         set_initial_concentrations(handle, *f_init(n_compartment))
#     else:
#         if not serde_path:
#             raise ValueError("serde_path must be provided when is_serde is True")
#         handle_module.register_serde(handle, serde_path)

#     # Apply and execute
#     ok, msg = handle_module.apply(handle, is_serde)
#     if not ok:
#         print(msg)
#         return -1

#     handle_module.exec(handle)
#     return handle_module.i_rank(handle)
