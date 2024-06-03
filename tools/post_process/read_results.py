from dataclasses import dataclass
import numpy as np
import h5py

@dataclass
class RawResults:
    distribution: np.ndarray
    initial_distribution: np.ndarray
    cliq: np.ndarray
    data: np.ndarray
    tf: np.ndarray
    dt: np.ndarray
    npart: np.ndarray
    records_distribution: np.ndarray
    t: np.ndarray
    n_t: int

def import_results(path):
    results = RawResults(distribution=None, initial_distribution=None, cliq=None, data=None, tf=None, dt=None,
                         npart=None, records_distribution=None, t=None, n_t=None)

    with h5py.File(path, 'r') as file:
        # Extract data from the HDF5 file
        results.distribution = np.array(file.get("final_results/distribution"))
        results.initial_distribution = np.array(file.get("initial_parameters/particle_distribution"))
        results.cliq = np.array(file.get("final_results/concentrations/liquid"))
        results.data = np.array(file.get("final_results/concentrations/records"))
        results.tf = np.array(file.get('initial_parameters/final_time'))
        results.dt = np.array(file.get('initial_parameters/delta_time'))
        results.npart = np.array(file['initial_parameters/number_particles_0'])
        results.records_distribution = np.array(file.get("final_results/concentrations/records_distribution"))
        results.t =  np.array(file.get("final_results/time"))
    

    results.t = results.t.reshape(-1,)
    
    # Calculate t and n_t
    n_t = results.records_distribution.shape[0]
    # results.t = np.linspace(0,results.tf,n_t)#np.array([i * results.dt for i in range(n_t)])
    results.n_t = n_t

    return results
