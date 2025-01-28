from .handle_module import *
import os

__all__ = []
__doc__ = handle_module.__doc__
if hasattr(handle_module, "__all__"):
    __all__ = handle_module.__all__


def pyinit_handle(sim_id):
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
    return handle_module.init_handle(n_rank, i_rank, sim_id, 1),i_rank,n_rank



__all__.extend([
    "pyinit_handle"
])