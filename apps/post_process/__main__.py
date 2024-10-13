from typing import List
from matplotlib import pyplot as plt
import numpy as np

from .rtd import get_scalar_rtd

# from .old_read_results import import_results
from .read_results import import_results, Results
from . import (
    FIGURE_TYPE,
    average_concentration,
    check_mixing,
    set_time_unit_to_hour,
    time_average_reaction_rate,
)
from .properties import process_particle_data, property_distribution
import os
import argparse
from numpy.polynomial import Polynomial
from . import get_time


def assemble(res_folder: str, names: List[str]) -> List[str]:
    return [f"{res_folder}{i}/{i}.h5" for i in names]


def mk_parser():
    parser = argparse.ArgumentParser(description="Post process tool")

    parser.add_argument(
        dest="name_results",
        type=str,
        nargs="+",
        help="List of names to generate result paths for.",
    )

    parser.add_argument(
        "--root_res",
        type=str,
        default="./results/",
        help="Root directory for the results.",
        required=False,
    )
    return parser


def detect_exponential_growth(t: np.ndarray, n: np.ndarray, threshold: float = 0.01):
    """
    Detect the index where exponential growth starts by analyzing the rate of change
    in the logarithmic domain.

    Parameters:
    t (np.ndarray): Time array
    n (np.ndarray): Particle number array
    threshold (float): Threshold for identifying the start of exponential growth

    Returns:
    int: Index where exponential growth starts
    """
    log_n = np.log(n)
    # Calculate the derivative (rate of change) of log(n) with respect to time
    rate_of_change = np.gradient(log_n, t)

    # Detect where the rate of change stabilizes (exponential region)
    for i in range(1, len(rate_of_change)):
        if np.abs(rate_of_change[i] - rate_of_change[i - 1]) < threshold:
            return i
    return 0  # Default to 0 if no stable growth is detected


def plot_grow_in_number(t: np.ndarray, n: np.ndarray, dest: str):
    n = np.sum(n, axis=1)
    plt.figure()
    plt.semilogy(t, n, "-", color="black", label="results")
    plt.ylabel("Number of particles")
    plt.xlabel(f"Time [{get_time()}]")
    plt.title("Particle number growth (Log-Scale)")
    plt.grid(True)

    # start_index = detect_exponential_growth(t, n,0.05)

    # exp_t = t[start_index:]
    # exp_n = n[start_index:]

    # fitting = Polynomial.fit(exp_t, np.log(exp_n), 1)

    # coeff = fitting.convert().coef

    # y_exp = np.exp(coeff[0]) * np.exp(coeff[1] * t)
    # plt.semilogy(t, y_exp, label="regression")
    # plt.legend()

    # Save the plot
    plt.savefig(f"{dest}/num_grow{FIGURE_TYPE}")


def check_time_unit(results: Results):
    # Conversion to hour if duration too long
    if results.main.time[-1] > 10000:
        results.main.time = results.main.time / 3600.0
        set_time_unit_to_hour()


def calculate_mass(dest:str,results: Results):
    if "mass" in results.partial[0].extra_bioparam[0].keys():
        cx = np.zeros((results.main.n_export,))
        for i in range(results.main.n_export):
            total_mass = 0
            for dataset in results.partial:
                mass = dataset.extra_bioparam[i]["mass"]
                total_mass += np.sum(mass)
            cx[i] = results.main.weight*total_mass/results.main.volume_liquid[i,0]
        
        print("X/X0: ",cx[-1],cx[0],cx[-1]/cx[0])
        plt.figure()
        plt.plot(results.time,cx)
        plt.savefig(dest+"/cx")
    pass


def get_spatial_average_concentration(results: Results, dest: str):
    def phase_functor(concentration_record, full_volume, phase_name: str):
        c_avg = np.array(average_concentration(concentration_record, full_volume))
        plt.figure()
        plt.plot(results.time, c_avg)
        plt.title("Average glucose concentration over time")
        plt.ylabel("C [g/L]")
        plt.xlabel(f"Time [{get_time()}]")
        plt.savefig(dest + f"/c_avg_{phase_name}_{FIGURE_TYPE}")

    concentration_record = results.main.concentrations_liquid[:, :, 0]
    full_volume = results.main.volume_liquid
    phase_functor(concentration_record, full_volume, "liquid")

    # if results.main.concentrations_gas is not None:
    #     concentration_record = results.main.concentrations_gas[:, :, 0]
    #     full_volume = results.main.volumes_gas
    #     phase_functor(concentration_record, full_volume, "gas")


def get_vtk(id: int):
    return None  # f"./results/{name_results[i]}/{name_results[i]}_{last_id}.vtu"


def time_average_reaction_data(results):
    # np.set_printoptions(precision=10)
    # print(
    #     time_average_reaction_rate(
    #         results.time[-1], results.data[:, :, 0], results.volume_liquid
    #     )
    # )
    pass


def main(name_results, root_res="./results/"):
    dest = [f"./results/{i}/postprocessing" for i in name_results]

    for d in dest:
        if not os.path.exists(d):
            os.makedirs(d)

    pathres = assemble(root_res, name_results)

    vtu_path = None  # ("/mnt/c/Users/casale/Documents/code/cpp/biomc/cma_data/bench/cma_mesh.vtu")
    check_mixing(name_results, pathres, dest, vtu_path)

    for i, current_path_res in enumerate(pathres):
        results: Results = import_results(current_path_res)

        check_time_unit(results)

        last_id = results.main.n_export - 1
        last_vtk_path = get_vtk(last_id)

        get_scalar_rtd(dest[i],0.5e-6,results,False) 
        # get_scalar_rtd(dest[i],0.035,results)
        # get_scalar_rtd(dest[i],0.031653119013143756,results,True)
        # get_scalar_rtd(dest[i],2.6e-5,results,True) 

        calculate_mass(dest[i],results)
        plot_grow_in_number(results.time, results.total_repartion, dest[i])
        get_spatial_average_concentration(results, dest[i])

        

        process_particle_data(results, dest[i])


if __name__ == "__main__":
    args = mk_parser().parse_args()

    name_results = args.name_results
    root_res = args.root_res
    main(name_results, root_res)
