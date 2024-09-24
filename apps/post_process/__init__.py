from typing import List
from matplotlib import pyplot as plt
import numpy as np
from typing import Optional, Tuple
import cmtool.vtk
from .read_results import RawResults, import_results
import os
from .initialiser import make_initial_concentration

RATIO_MASS_LENGTH = 0.45044876111074444 / 0.12477411510047276

FIGURE_TYPE=".png"
TIME_UNIT ="s"

def set_time_unit_to_hour():
    global TIME_UNIT
    TIME_UNIT="h"

def get_time():
    return TIME_UNIT

def mkdir(d):
    if not os.path.exists(d):
        os.makedirs(d)


def check_mixing(
    name_results, pathres: List[str], dest: str, vtk_cma_mesh_path: Optional[str] = None
):
    for i in range(len(pathres)):
        results = import_results(pathres[i])
        if results is None:
            break
        (
            normalized_scalar_concentration,
            norm_c_var,
            normalized_particle_concentration,
            norm_par_var,
            t,
        ) = process_norm(name_results[i], results, vtk_cma_mesh_path)
        plt.scatter(t, norm_par_var, label=name_results[i])
        plt.plot(t, norm_c_var, label=f"liquid_{name_results[i]}")

    plt.legend()
    plt.title("Segregation index as a function of the time")
    plt.ylabel(r"\[ \frac{\sigma(t)}{\sigma(t_{0})}\]")
    plt.xlabel(f"time [{get_time()}]")
    for i in dest:
        plt.savefig(f"{i}/mixing_variance{FIGURE_TYPE}", dpi=1500)


def norm_concentration(
    raw_concentration: np.ndarray, volumes: np.ndarray
) -> Tuple[np.ndarray, float, float]:
    vtot = np.sum(volumes, axis=1)
    mean_concentration = np.sum(raw_concentration * volumes, axis=1) / vtot
    mean_concentration = mean_concentration.reshape(-1, 1)
    variance = (
        np.sum(np.power(raw_concentration - mean_concentration, 2) * volumes, axis=1)
        / vtot
    )
    return raw_concentration / mean_concentration, mean_concentration, variance


def process_norm(
    name: str, results: RawResults, vtk_cma_mesh_path: Optional[str] = None
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    concentration_record = results.data[:, :, 0]

    full_volume = results.volume_liquid

    p_concentration = results.distribution / full_volume

    normalized_scalar_concentration, _, var_c = norm_concentration(
        concentration_record, full_volume
    )
    normalized_particle_concentration, _, par_var = norm_concentration(
        p_concentration, full_volume
    )

    norm_par_var = par_var / par_var[0]
    norm_c_var = var_c / var_c[0]

    if vtk_cma_mesh_path is not None:
        cmtool.vtk.mk_series(
            vtk_cma_mesh_path,
            "./results/",
            name,
            results.t,
            [full_volume, "liquid_volume"],
            [p_concentration, "particle_concentration"],
            [normalized_particle_concentration, "normalized_particle_concentration"],
            [normalized_scalar_concentration, "normalized_liquid_concentration"],
            [concentration_record, "liquid_concentration"],
        )

    return (
        normalized_scalar_concentration,
        norm_c_var,
        normalized_particle_concentration,
        norm_par_var,
        results.t,
    )


def average_concentration(results: RawResults):
    concentration_record = results.data[:, :, 0]

    full_volume = results.volume_liquid

    c_avg = np.sum(concentration_record * full_volume, axis=1) / np.sum(
        full_volume, axis=1
    )
    return c_avg


def time_average_reaction_rate(
    duration: float, time_evolution_data: np.ndarray, time_evolution_volume
):
    def mass_func(i, x, v):
        return np.sum(x[i, :] * v[i, :])

    # gram version x1000
    m_init = mass_func(0, time_evolution_data, time_evolution_volume) * 1e3
    m_end = mass_func(-1, time_evolution_data, time_evolution_volume) * 1e3
    return (m_end - m_init) / duration, m_init, m_end
