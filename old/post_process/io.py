import cmtool.vtk
from typing import List
import numpy as np

def append_resukts_scalar_vtk(filename: str, value: np.ndarray, name: str):
    scalar = cmtool.vtk.mk_scalar(value, name)
    cmtool.vtk.append_scalar(filename, filename, scalar)
    pass

def write_vtk_results(
    vtk_cma_mesh_path: str,
    sname: List[str],
    t: float,
    full_volume: np.ndarray,
    p_concentration: np.ndarray,
    normalized_particle_concentration: np.ndarray,
    normalized_scalar_concentration: np.ndarray,
    concentration_record: np.ndarray,
):
    cmtool.vtk.mk_series(
        vtk_cma_mesh_path,
        "./results/",
        sname,
        t,
        [full_volume, "liquid_volume"],
        [p_concentration, "particle_concentration"],
        [normalized_particle_concentration, "normalized_particle_concentration"],
        [normalized_scalar_concentration, "normalized_liquid_concentration"],
        [concentration_record, "liquid_concentration"],
    )
