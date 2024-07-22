from typing import Dict, List, Optional
from matplotlib import pyplot as plt
import numpy as np
from apps.post_process.read_results import import_results
from apps.post_process import average_concentration, process_norm
from .properties import mk_histogram, get_distribution_moment
import os
import cmtool.vtk


def check_mixing(
    name_results, pathres: List[str], dest: str, vtk_cma_mesh_path: Optional[str] = None
):
    for i in range(len(pathres)):
        results = import_results(pathres[i])
        (
            normalized_scalar_concentration,
            norm_c_var,
            normalized_particle_concentration,
            norm_par_var,
            t,
        ) = process_norm(name_results[i], results, vtk_cma_mesh_path)
        plt.scatter(t, norm_par_var, label=name_results[i])
        plt.semilogy(t, norm_c_var, label=f"liquid_{name_results[i]}")

    plt.legend()
    plt.title("Segregation index as a function of the time")
    plt.ylabel(r"\[ \frac{\sigma(t)}{\sigma(t_{0})}\]")
    plt.xlabel("time [s]")
    plt.savefig(f"{dest}/mixing_variance.svg", dpi=1500)


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
            mass = np.sum(value)/len(value)
            print("mass: ", mass)

        if isinstance(value, np.ndarray) and np.issubdtype(value.dtype, float):
            mean, variance_population, variance_sample = get_distribution_moment(value)

            print(mean, variance_population, variance_sample)

            mk_histogram(value, f"{prefix}_{key}", dest)
            if vtk_cma_mesh_path is not None:
                append_resukts_scalar_vtk(vtk_cma_mesh_path, value, key)


def assemble(res_folder: str, names: List[str]) -> List[str]:
    return [f"{res_folder}{i}.h5" for i in names]


if __name__ == "__main__":
    root_res = "./results/"
    dest = "./results/monod_0d_/pp"

    if not os.path.exists(dest):
        os.makedirs(dest)

    name_results = ["monod_0d"]
    pathres = assemble(root_res, name_results)

    vtu_path = None  # ("/mnt/c/Users/casale/Documents/code/cpp/biomc/cma_data/bench/cma_mesh.vtu")
    check_mixing(name_results, pathres, dest, vtu_path)

    for i, p in enumerate(pathres):
        results = import_results(p)
        last_id = results.n_t - 1
        last_vtk_path = (
            None  # f"./results/{name_results[i]}/{name_results[i]}_{last_id}.vtu"
        )
        property_distribution(results.initial_bioparam, f"{name_results[i]}_init", dest)
        property_distribution(
            results.final_bioparam, f"{name_results[i]}_final", dest, last_vtk_path
        )
        np.set_printoptions(precision=10)
        plt.figure()
        c_avg = np.array(average_concentration(results))
        plt.plot(results.t, c_avg)
        plt.savefig(dest + "/c_avg.png")
