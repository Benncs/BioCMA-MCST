#!/usr/bin/python3

import argparse
import os
import sys
import subprocess
import re
import time

from cli_formater import format_cli


__current_file_path = os.path.abspath(__file__)
__current_directory = os.path.dirname(__current_file_path)
ROOT = __current_directory + "/.."
DEFAULT_TYPE = "debugoptimized"
MPI_COMMAND = "mpiexec --allow-run-as-root -np 4 "
OMP_NUM_THREADS = "1"


def get_executable(type: str):
    return f"{ROOT}/builddir/apps/cli/biocma_mcst_cli_app"

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


def format_rhs(match):
    """
    Function to format RHS values in green.

    Args:
        match (re.Match): Match object containing the regex match.

    Returns:
        str: Formatted string with ANSI escape codes for green text.
    """
    return f"{match.group(1)}\033[92m{match.group(2)}\033[0m"


def exec(command, n_thread):
    env_var = os.environ.copy()
    env_var["OMP_PLACES"] = "threads"
    env_var["OMP_PROC_BIND"] = "spread"
    env_var["OMP_NUM_THREADS"] = n_thread

    result = command.replace("-", "\n-")
    pattern = re.compile(r"(-\w+\s)(\S+)")
    formatted_command = pattern.sub(format_rhs, result)
    print("\r\n")
    print(formatted_command)
    print("\n")
    start_time = time.perf_counter()
    process = subprocess.Popen(command, shell=True, env=env_var)
    return_code = process.wait()
    end_time = time.perf_counter()
    elapsed_time = end_time - start_time
    print(f"Command executed in \033[92m{elapsed_time:.6f}\033[0m seconds")
    return return_code


def main():

    cli_args = mk_parser().parse_args()
    r_type = DEFAULT_TYPE
    if cli_args.release_flag:
        r_type = "release"

    run_cli = format_cli(["_", cli_args.casename])

    mpi_c = ""
    if cli_args.use_mpi:
        mpi_c = MPI_COMMAND + " "

    command = (
        mpi_c + get_executable(r_type) + " " + run_cli + f" -nt {cli_args.n_threads}" + '--device_id=1'
    )
    exec(command, cli_args.n_threads)

if __name__ == "__main__":
    main()




