from dataclasses import dataclass
import numpy as np
import h5py
from typing import Dict, Optional, List


@dataclass
class RawResults:
    """
    Dataclass to hold the results imported from an HDF5 file.

    Attributes:
    distribution (Optional[np.ndarray]): Final particle distribution.
    initial_distribution (Optional[np.ndarray]): Initial particle distribution.
    cliq (Optional[np.ndarray]): Liquid concentration.
    data (Optional[np.ndarray]): Records of concentrations.
    tf (Optional[np.ndarray]): Final time.
    dt (Optional[np.ndarray]): Delta time.
    npart (Optional[np.ndarray]): Initial number of particles.
    records_distribution (Optional[np.ndarray]): Records of distribution.
    t (Optional[np.ndarray]): Time array.
    n_t (Optional[int]): Number of time records.
    """

    version: int

    distribution: Optional[np.ndarray] = None
    initial_distribution: Optional[np.ndarray] = None
    cliq: Optional[np.ndarray] = None
    data: Optional[np.ndarray] = None
    tf: Optional[np.ndarray] = None
    dt: Optional[np.ndarray] = None
    npart: Optional[np.ndarray] = None
    volume_liquid: Optional[np.ndarray] = None
    volume_gas: Optional[np.ndarray] = None
    records_distribution: Optional[np.ndarray] = None
    t: Optional[np.ndarray] = None
    n_t: Optional[int] = None
    final_bioparam: Optional[Dict[str, np.ndarray]] = None
    initial_bioparam: Optional[Dict[str, np.ndarray]] = None
    t_per_flow_map: Optional[float] = None
    n_per_flow_map: Optional[int] = None
    extra_bioparam: Optional[List[Dict[str, np.ndarray]]] = None
    weight: Optional[float] = None


def __import_v1(file):
    results = RawResults(
        version=1,
        distribution=None,
        initial_distribution=None,
        cliq=None,
        data=None,
        tf=None,
        dt=None,
        npart=None,
        records_distribution=None,
        t=None,
        n_t=None,
    )
    # Extract data from the HDF5 file
    results.distribution = np.array(file.get("final_results/distribution"))
    results.initial_distribution = np.array(
        file.get("initial_parameters/particle_distribution")
    )
    results.cliq = np.array(file.get("final_results/concentrations/liquid"))
    results.data = np.array(file.get("final_results/concentrations/records"))
    results.tf = np.array(file.get("initial_parameters/final_time"))
    results.dt = np.array(file.get("initial_parameters/delta_time"))
    results.npart = np.array(file["initial_parameters/number_particles_0"])
    results.records_distribution = np.array(
        file.get("final_results/concentrations/records_distribution")
    )
    results.t = np.array(file.get("final_results/time"))

    # Ensure t is a 1D array
    if results.t is not None:
        results.t = results.t.reshape(-1)

    # Calculate n_t if records_distribution is available
    if results.records_distribution is not None:
        results.n_t = results.records_distribution.shape[0]

    return results


def __import_v2(file):
    results = RawResults(
        version=2,
        distribution=None,
        initial_distribution=None,
        cliq=None,
        data=None,
        tf=None,
        dt=None,
        npart=None,
        records_distribution=None,
        t=None,
        n_t=None,
    )
    results.distribution = np.array(file.get("records/distribution"))
    results.volume_gas = np.array(file.get("records/gas_volume"))
    results.volume_liquid = np.array(file.get("records/liquid_volume"))
    results.data = np.array(file.get("records/concentration_liquid"))
    results.tf = np.array(file.get("initial_parameters/final_time"))
    results.dt = np.array(file.get("initial_parameters/delta_time"))
    results.records_distribution = np.array(file.get("records/distribution"))
    results.t = np.array(file.get("records/time"))

    # Calculate n_t if records_distribution is available
    if results.records_distribution is not None:
        results.n_t = results.records_distribution.shape[0]
    return results


def __import_v3(file):
    results = RawResults(
        version=3,
        distribution=None,
        initial_distribution=None,
        cliq=None,
        data=None,
        tf=None,
        dt=None,
        npart=None,
        records_distribution=None,
        t=None,
        n_t=None,
    )
    results.distribution = np.array(file.get("records/number_particles"))
    results.volume_gas = np.array(file.get("records/gas_volume"))
    results.volume_liquid = np.array(file.get("records/liquid_volume"))
    results.data = np.array(file.get("records/concentration_liquid"))
    results.tf = np.array(file.get("initial_parameters/final_time"))
    results.dt = np.array(file.get("initial_parameters/delta_time"))
    results.records_distribution = np.array(file.get("records/number_particles"))
    results.t = np.array(file.get("records/time"))
    results.t_per_flow_map = np.array(file.get("initial_parameters/t_per_flow_map"))
    final_bio = file.get("biological_model/final", None)
    init_bio = file.get("biological_model/initial", None)
    results.npart = np.array(file.get("final_results/number_particles"))

    results.weight = np.array(file.get("initial_parameters/initial_weight"))

    i_user_export = 1
    user_export = []

    while True:
        node = file.get(f"biological_model/{i_user_export}", None)
        if node is None:
            break
        d = {}

        for key in node:
            d[key] = np.array(node[key])

        user_export.append(d)
        i_user_export += 1
    results.extra_bioparam = user_export

    if results.dt is not None and results.t_per_flow_map is not None:
        results.n_per_flow_map = int(results.t_per_flow_map / results.dt)

    if final_bio is not None:
        results.final_bioparam = {}
        for key in final_bio:
            results.final_bioparam[key] = np.array(final_bio[key])

    if init_bio is not None:
        results.initial_bioparam = {}
        for key in init_bio:
            results.initial_bioparam[key] = np.array(init_bio[key])

    # Calculate n_t if records_distribution is available
    if results.records_distribution is not None:
        results.n_t = results.records_distribution.shape[0]
    return results


def import_results(path: str) -> RawResults:
    """
    Import results from an HDF5 file and store them in a RawResults dataclass.

    Parameters:
    path (str): Path to the HDF5 file.

    Returns:
    RawResults: Dataclass containing the imported results.

    Raises:
    FileNotFoundError: If the file at the specified path is not found.
    OSError: If there is an error opening the file.
    Exception: For any other unexpected errors.
    """

    try:
        with h5py.File(path, "r") as file:
            _attributes = list(file.attrs)

            version = file.attrs.get("file_version", -1)
            if version == 1:
                return __import_v1(file)
            elif version == 2:
                return __import_v2(file)
            elif version == 3:
                return __import_v3(file)
            else:
                raise Exception("Unknown file version")

    except FileNotFoundError:
        raise FileNotFoundError(f"The file at path {path} was not found.")
    except OSError:
        raise OSError(
            f"Error opening the file at path {path}. It might be corrupted or in use by another process."
        )
    except Exception as e:
        raise Exception(f"An unexpected error occurred: {e}")
