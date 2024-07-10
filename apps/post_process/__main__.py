from typing import Dict, List, Optional
from matplotlib import pyplot as plt
import numpy as np
from apps.post_process.read_results import RawResults, import_results
from apps.post_process import process_norm
from .properties import mk_histogram,get_distribution_moment

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


def property_distribution(biodict:Dict[str, np.ndarray],prefix:str=""):
    for key in biodict:
    
        value = biodict[key]
        
        if isinstance(value, np.ndarray) and np.issubdtype(value.dtype, float):
        
            mean,variance_population,variance_sample = get_distribution_moment(value)

            print(mean,variance_population,variance_sample)

            mk_histogram(value,f"{prefix}_{key}")

    


def assemble(res_folder: str, names: List[str]) -> List[str]:
    return [f"{res_folder}{i}.h5" for i in names]


if __name__ == "__main__":
    root_res = "./results/"
    name_results = ["test_model_fed_out"]
    pathres = assemble(root_res, name_results)

    vtu_path = (
        "/mnt/c/Users/casale/Documents/code/cpp/biomc/cma_data/bench/cma_mesh.vtu"
    )
    check_mixing(name_results, pathres, root_res,vtu_path)
    results = import_results(pathres[0])
    property_distribution(results.initial_bioparam,"init")
    property_distribution(results.final_bioparam,"final")