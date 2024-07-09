import os
import sys 
import subprocess
import re
import time

from cli_formater import format_cli
def get_executable(type:str):
    return f"./builddir/{type}/apps/cli/biocma_mcst_cli_app"

DEFAULT_TYPE="debugoptmized"
MPI_COMMAND="mpiexec --allow-run-as-root "

OMP_NUM_THREADS = 12

def format_rhs(match):
    """
    Function to format RHS values in green.

    Args:
        match (re.Match): Match object containing the regex match.

    Returns:
        str: Formatted string with ANSI escape codes for green text.
    """
    return f"{match.group(1)}\033[92m{match.group(2)}\033[0m"

def parse_cli(args):
    cli_args = {"type": DEFAULT_TYPE,"name":None}

    if len(args) < 2:
        print("Usage: runner.py <type> [-r]")
        exit(0)

    cli_args["name"] = args[1]

    if len(args) == 3 and args[2] == '-r':
        cli_args["type"] = "release"

    return cli_args


def exec(command):
    env_var = os.environ.copy()
    env_var["OMP_NUM_THREADS"] = str(OMP_NUM_THREADS)


    result = command.replace("-", "\n-")
    pattern = re.compile(r'(-\w+\s)(\S+)')
    formatted_command = pattern.sub(format_rhs, result)
    print(formatted_command)
    print("\n")
    start_time = time.perf_counter()
    process = subprocess.Popen(command, shell=True,env=env_var)
    return_code = process.wait()
    end_time = time.perf_counter()
    elapsed_time = end_time - start_time
    print(f"Command executed in \033[92m{elapsed_time:.6f}\033[0m seconds")
    return return_code

if __name__=="__main__":
    args = sys.argv 
    cli_args = parse_cli(args)

    run_cli = format_cli(["_",cli_args["name"]])
    use_mpi=False
    mpi_c = ""
    if use_mpi:
        mpi_c = MPI_COMMAND

    command = mpi_c+get_executable(cli_args["type"])+" "+run_cli

    exec(command)

    

    
