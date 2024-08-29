from typing import Dict, List, Optional
from matplotlib import pyplot as plt
import numpy as np
from apps.post_process.read_results import import_results
from apps.post_process import (
    average_concentration,
    check_mixing,
    process_norm,
    time_average_reaction_rate,
)
from .properties import mk_histogram, get_distribution_moment
import os
import cmtool.vtk
import argparse

RATIO_MASS_LENGTH = 0.45044876111074444 / 0.12477411510047276


def append_resukts_scalar_vtk(filename: str, value: np.ndarray, name: str):
    scalar = cmtool.vtk.mk_scalar(value, name)
    cmtool.vtk.append_scalar(filename, filename, scalar)
    pass


def property_distribution(
    biodict: Dict[str, np.ndarray],
    prefix: str = "",
    dest: str = "./results/",
    vtk_cma_mesh_path: Optional[str] = None,
):
    for key in biodict:
        value = biodict[key]

        if key == "lenght":
            mass = np.sum(value) * RATIO_MASS_LENGTH
            print("mass: ", mass)

        if isinstance(value, np.ndarray) and np.issubdtype(value.dtype, float):
            mean, variance_population, variance_sample = get_distribution_moment(value)

            print(
                key,
                ": ",
                "mean: ",
                mean,
                "var: ",
                variance_population,
                "varred: ",
                variance_sample,
            )

            mk_histogram(value, f"{prefix}_{key}", dest)
            if vtk_cma_mesh_path is not None:
                append_resukts_scalar_vtk(vtk_cma_mesh_path, value, key)


def property_space(i: int, biodict: Dict[str, np.ndarray], key1: str, key2: str):
    value1 = biodict[key1]
    value2 = biodict[key2]
    MAX_SAMPLE = 1_000_000
    sample_size = min(len(value1), MAX_SAMPLE)  # or any smaller number
    idx = np.random.choice(range(len(value1)), size=sample_size, replace=False)

    if isinstance(value1, np.ndarray) and np.issubdtype(value1.dtype, float):
        if isinstance(value2, np.ndarray) and np.issubdtype(value2.dtype, float):
            plt.scatter(value1[idx], value2[idx], label=f"data {i}", s=1)


def plot_property_space(
    biodicts: List[Dict[str, np.ndarray]],
    key1: str,
    key2: str,
    prefix: str = "",
    dest: str = "./results/",
):
    for i, d in enumerate(biodicts):
        plt.figure()
        property_space(i, d, key1, key2)
        plt.xlabel(key1)
        plt.ylabel(key2)
        plt.legend()
        plt.savefig(f"{dest}/plot_{key1}_{key2}_{i}")


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


if __name__ == "__main__":
    args = mk_parser().parse_args()

    name_results = args.name_results
    root_res = args.root_res

    dest = [f"./results/{i}/pp" for i in name_results]

    for d in dest:
        if not os.path.exists(d):
            os.makedirs(d)

    pathres = assemble(root_res, name_results)

    vtu_path = None  # ("/mnt/c/Users/casale/Documents/code/cpp/biomc/cma_data/bench/cma_mesh.vtu")
    check_mixing(name_results, pathres, dest, vtu_path)
    X0 = 0.1
    for i, p in enumerate(pathres):
        results = import_results(p)
        last_id = results.n_t - 1
        last_vtk_path = (
            None  # f"./results/{name_results[i]}/{name_results[i]}_{last_id}.vtu"
        )

        init_mass = np.sum(X0 * results.volume_liquid[0, :])
        print("INITIAL MASS", init_mass)
        final_mass = results.npart*results.weight/np.sum(results.volume_liquid[0, :])
        print("FINAL CONCENTRATION",final_mass)

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

        plot_property_space(
            dict_particles, "mu", "lenght", f"{name_results[i]}", dest[i]
        )

        np.set_printoptions(precision=10)
        plt.figure()
        c_avg = np.array(average_concentration(results))
        print(
            time_average_reaction_rate(
                results.t[-1], results.data[:, :, 0], results.volume_liquid
            )
        )
        plt.plot(results.t, c_avg)
        plt.savefig(dest[i] + "/c_avg.png")
