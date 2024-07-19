import numpy as np
from typing import Optional, Tuple
import cmtool.vtk
from .read_results import RawResults


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

    c_avg = np.sum(concentration_record * full_volume,axis=1) / np.sum(full_volume, axis=1)
    return c_avg
