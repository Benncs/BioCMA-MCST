import birem 
import numpy as np 
import birem.vtk 
from read_results import * 
import matplotlib.pyplot as plt 
import scienceplots
from typing import Optional,List
# plt.style.use(['science','ieee','scatter','grid'])
plt.style.use(['science','grid'])

def write_vtk_results(vtk_cma_mesh_path:str,sname:List[str],t:float,
                        full_volume:np.ndarray,
                        p_concentration:np.ndarray,
                        normalized_particle_concentration:np.ndarray,
                        normalized_scalar_concentration:np.ndarray,
                        concentration_record:np.ndarray):
    birem.vtk.mk_series(vtk_cma_mesh_path,"./results/",sname,t,
        [full_volume,"liquid_volume"],[p_concentration,"particle_concentration"],[normalized_particle_concentration,"normalized_particle_concentration"],[normalized_scalar_concentration,"normalized_liquid_concentration"],[concentration_record,"liquid_concentration"])


def norm_concentration(raw_concentration,volumes):
    
    vtot = np.sum(volumes,axis=1)
    mean_concentration = np.sum(raw_concentration * volumes, axis=1) / vtot
    mean_concentration = mean_concentration.reshape(-1, 1)
    variance = np.sum(np.power(raw_concentration - mean_concentration, 2) * volumes, axis=1) / vtot
    return raw_concentration / mean_concentration,mean_concentration,variance

def read_compute_norm(pathres,folder_root,sname,vtk_cma_mesh_path:Optional[str]=None):
    # Import results from the given path
    results = import_results(pathres)

    # Initialize array to store volume data
    read_volume = np.zeros((14, results.records_distribution.shape[1]))

    # Define the flow map time step and number of steps per flow map
    t_per_flowmap = 0.0286
    n_per_flowmap = int(t_per_flowmap / results.dt)

    # Read volume data from files and store in read_volume
    for i in range(14):
        folder = f"{folder_root}/bench_2/i_{i+1}"
        raw_volume = birem.read_scalar(f"{folder}/vofL.raw")
        read_volume[i] = raw_volume
    

    # Reshape the records array
    records_p_dist = results.records_distribution.reshape(results.records_distribution.shape[0], results.records_distribution.shape[1])

    # Calculate the total number of particles at the initial time step
    n_particle = np.sum(records_p_dist[0])

    # Initialize arrays for particle concentration and full volume
    p_concentration = np.zeros((records_p_dist.shape[0], records_p_dist.shape[1]))
    full_volume = np.zeros_like(p_concentration)

    # Set initial full volume and particle concentration
    full_volume[0] = np.copy(read_volume[0])
    p_concentration[0] = records_p_dist[0] / read_volume[0]

    # Initialize array for concentration record and mass
    concentration_record = results.data[:, :, 0]

    
    # Loop over each time step to compute volumes and concentrations
    for i in range(0, results.n_t):
        index = (i // n_per_flowmap) % 14
        full_volume[i] = np.copy(read_volume[index])
        p_concentration[i] = records_p_dist[i] / read_volume[index]

    vtot = np.sum(full_volume,axis=1)
   

    
    normalized_scalar_concentration,mean_c,var_c = norm_concentration(concentration_record,full_volume) 
    normalized_particle_concentration ,mean_p_c,par_var= norm_concentration(p_concentration,full_volume) 

    if vtk_cma_mesh_path is not None:
        write_vtk_results(vtk_cma_mesh_path,sname,results.t,full_volume,p_concentration,normalized_particle_concentration,normalized_scalar_concentration,concentration_record)
      
    return normalized_particle_concentration,normalized_scalar_concentration,par_var/par_var[0],var_c/var_c[0],results.t

folder_root = "/home/benjamin/Documenti/code/cpp/biomc/cma_data/"

sname = ["data_10M_n"]
pathres = ["./results/10M.h5"]
#pathres = ['./results/50k.h5',"./results/1M.h5","./results/5M.h5"]
# pathres = [ "/home/benjamin/Documenti/code/cpp/biomc/results/result_2024-06-04-11:14:49.h5"]
# sname = ["test_linterp"]
vtk_cma_mesh_path="/home/benjamin/Documenti/code/cpp/BIREM_new/out/sanofi/cma_mesh.vtu"


_,n_c ,_,pvc,t= read_compute_norm(pathres[0],folder_root,sname[0])
# plt.semilogy(t,pvc,label="liquid")
for i in range(len(pathres)):
    n,n_c,par_var,_,t = read_compute_norm(pathres[i],folder_root,sname[i],vtk_cma_mesh_path)
#     plt.semilogy(t, par_var,label=sname[i])




# plt.legend()
# plt.title("Segregation index as a function of the time")
# plt.ylabel(r"\[ \frac{\sigma(t)}{\sigma(t_{0})}\]")
# plt.xlabel("time [s]")
# # plt.savefig("./results/mixing_variance3.svg",dpi=1500)
# plt.show()


