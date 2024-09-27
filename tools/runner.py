#!/usr/bin/python3

import argparse
import os
from cli_formater import format_cli
import sys

dev_path = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "devutils"))
sys.path.append(dev_path)

from exec import exec    # noqa: E402


__current_file_path = os.path.abspath(__file__)
__current_directory = os.path.dirname(__current_file_path)
ROOT = __current_directory + "/.."
DEFAULT_TYPE = "debug"
MPI_COMMAND = "mpiexec --allow-run-as-root -np 2 "
OMP_NUM_THREADS = "1"
COMPILER_NAME="gcc"

def get_executable(type: str, mpi: bool = True):
    appname = "biocma_mcst_cli_app" if mpi else "biocma_mcst_cli_app_shared"
    return f"{ROOT}/builddir/{type}_{COMPILER_NAME}/apps/cli/{appname}"


def mk_parser():
    parser = argparse.ArgumentParser(description="Runner")

    # Positional argument for casename
    parser.add_argument(dest="casename", type=str, help="The name of the case")

    # Optional flag -r
    parser.add_argument(
        "-r",
        dest="release_flag",
        action="store_true",
        required=False,
        help="Optional release",
    )

    parser.add_argument(
        "-p",
        dest="post_process",
        action="store_true",
        required=False,
        help="Do post process after",
    )


    # Optional number argument
    parser.add_argument(
        "-n",
        dest="n_threads",
        type=str,
        required=False,
        default=OMP_NUM_THREADS,
        help="Optional number of thread",
    )

    # Optional flag -mpi
    parser.add_argument(
        "-mpi",
        dest="use_mpi",
        action="store_true",
        help="Optional flag for MPI",
        required=False,
    )

    return parser


def main():
    cli_args = mk_parser().parse_args()
    r_type = DEFAULT_TYPE
    if cli_args.release_flag:
        r_type = "release"

    run_cli,res_file = format_cli(["_", cli_args.casename])
    mpi_c = ""
    if cli_args.use_mpi:
        mpi_c = MPI_COMMAND + " "

    command = (
        mpi_c
        + get_executable(r_type, cli_args.use_mpi)
        + " "
        + run_cli
        + f" -nt {cli_args.n_threads}"
    )
    exec(command, cli_args.n_threads,do_kokkos_measure=True)

  


if __name__ == "__main__":
    main()
