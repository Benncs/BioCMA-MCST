import handle_module
import os
from mpi4py import MPI

def init():
    os.environ["OMP_PROC_BIND"] = "spread"
    os.environ["OMP_PLACE"] = "threads"
    comm = MPI.COMM_WORLD
    n_rank =comm.Get_size()
    i_rank =comm.Get_rank()

    return handle_module.init_handle_raw(n_rank,i_rank)
