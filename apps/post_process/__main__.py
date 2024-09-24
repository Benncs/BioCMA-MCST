from typing import List
from matplotlib import pyplot as plt
import numpy as np
from .read_results import import_results
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
    return [f"{res_folder}{i}.h5" for i in names]


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
        if np.abs(rate_of_change[i] - rate_of_change[i-1]) < threshold:
            return i
    return 0  # Default to 0 if no stable growth is detected

def plot_grow_in_number(t: np.ndarray, n: np.ndarray, dest: str):





    n = np.sum(n,axis=1)
   
    plt.figure()
    plt.semilogy(t, n, "-",color='black', label="results")
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
    plt.show()

def main(name_results,root_res='./results/'):
    dest = [f"./results/{i}/postprocessing" for i in name_results]

    for d in dest:
        if not os.path.exists(d):
            os.makedirs(d)

    pathres = assemble(root_res, name_results)

    vtu_path = None  # ("/mnt/c/Users/casale/Documents/code/cpp/biomc/cma_data/bench/cma_mesh.vtu")
    check_mixing(name_results, pathres, dest, vtu_path)
    X0 = 0.1
    for i, current_path_res in enumerate(pathres):
        results = import_results(current_path_res)

        #Conversion to hour 
        
        if results.t[-1]>10000:
             results.t= results.t/3600
             set_time_unit_to_hour()

        last_id = results.n_t - 1
        last_vtk_path = (
            None  # f"./results/{name_results[i]}/{name_results[i]}_{last_id}.vtu"
        )
        try:
            init_mass = np.sum(X0 * results.volume_liquid[0, :])
            print("INITIAL MASS", init_mass)
            final_mass = (
                results.npart * results.weight / np.sum(results.volume_liquid[0, :])
            )
            print("FINAL CONCENTRATION", final_mass)
        except:
            pass
        property_distribution(
            results.initial_bioparam, f"{name_results[i]}_init", dest[i]
        )

        if results.final_bioparam is not None:
            property_distribution(
                results.final_bioparam,
                f"{name_results[i]}_final",
                dest[i],
                last_vtk_path,
            )

        dict_particles = [results.initial_bioparam]
        if results.extra_bioparam is not None:
            dict_particles = [*dict_particles, *results.extra_bioparam]
        if results.final_bioparam is not None:
            dict_particles.append(results.final_bioparam)

        process_particle_data(results.t,dict_particles, dest[i])
        plot_grow_in_number(results.t,results.records_distribution,dest[i])
        np.set_printoptions(precision=10)
        plt.figure()
        c_avg = np.array(average_concentration(results))
        print(
            time_average_reaction_rate(
                results.t[-1], results.data[:, :, 0], results.volume_liquid
            )
        )
        plt.plot(results.t, c_avg)
        plt.title("Average glucose concentration over time")
        plt.ylabel("C [g/L]")
        plt.xlabel(f"Time [{get_time()}]")
        plt.savefig(dest[i] + f"/c_avg{FIGURE_TYPE}")


if __name__ == "__main__":
    args = mk_parser().parse_args()

    name_results = args.name_results
    root_res = args.root_res
    main(name_results,root_res)
    
