#!/usr/bin/python3

import argparse
import os
import sys


def detect_execution_mode():
    release_app = "/opt/biomc/biocma_mcst_cli_app"
    if os.path.exists(release_app):
        return True
    return False


INSTALL_RELEASE = detect_execution_mode()

if INSTALL_RELEASE:
    gi = os.path.abspath(os.path.join("/opt/biomc"))
    sys.path.append(gi)  # TODO IMPROVE IMPORT
else:
    dev_path = os.path.abspath(os.path.join(os.path.dirname(__file__)))
    sys.path.append(dev_path)


from exec import exec  # noqa: E402
from cli_formater import format_cli  # noqa: E402

__current_file_path = os.path.abspath(__file__)
__current_directory = os.path.dirname(__current_file_path)
ROOT = __current_directory + "/.."
_MPI_ROOT_FLAG = ""  # "--allow-run-as-root"
MPI_COMMAND = f"mpiexec {_MPI_ROOT_FLAG} --mca orte_base_help_aggregate 1 -np 6  --report-bindings --bind-to core --map-by slot:PE=1"
# MPI_COMMAND = f"mpiexec {_MPI_ROOT_FLAG} -np 8 --use-hwthread-cpus"
OMP_NUM_THREADS = "1"


def get_executable(instal: str, mpi: bool = True):
    """Retourne le chemin de l'exécutable en fonction du type (debug/release)"""
    appname = "biocma_mcst_cli_app" if mpi else "biocma_mcst_cli_app_shared"

    if instal:
        return f"/opt/biomc/{appname}"
    else:
        return f"{ROOT}/builddir/host/apps/cli/{appname}"


def mk_parser():
    parser = argparse.ArgumentParser(description="Runner")

    parser.add_argument(dest="casename", type=str, help="The name of the case")

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

    parser.add_argument(
        "-n",
        dest="n_threads",
        type=str,
        required=False,
        default=OMP_NUM_THREADS,
        help="Optional number of thread",
    )

    parser.add_argument(
        "-mpi",
        dest="use_mpi",
        action="store_true",
        help="Optional flag for MPI",
        required=False,
    )

    parser.add_argument(
        "-s",
        dest="serde",
        action="store_true",
        help="serde",
        required=False,
    )
    parser.add_argument(
        "--dry-run",
        dest="dry_run",
        action="store_true",
        help="dry_run",
        required=False,
    )
    return parser


def main():
    cli_args = mk_parser().parse_args()

    run_cli, res_file = format_cli(
        ["_", cli_args.casename], global_install=INSTALL_RELEASE
    )

    mpi_c = ""
    if cli_args.use_mpi and not cli_args.dry_run:
        mpi_c = MPI_COMMAND + " "

    command = (
        mpi_c
        + get_executable(INSTALL_RELEASE, cli_args.use_mpi)
        + " "
        + run_cli
        + f"-nt {cli_args.n_threads} "
        + "-force 1 "
    )

    if cli_args.serde is True:
        command += "-serde " + "./results/prepbatch/prepbatch_serde_ "

        # if(cli_args.use_mpi):
        #     input("confirm force?")
        #     command+= " -force 1"
    if cli_args.dry_run:
        arg = command
        print(arg)
        return

    exec(command, cli_args.n_threads, do_kokkos_measure=False)


if __name__ == "__main__":
    main()
