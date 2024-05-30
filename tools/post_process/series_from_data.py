import birem 
import numpy as np 
import birem.vtk 
from read_results import * 
import matplotlib.pyplot as plt 
import scienceplots
# plt.style.use(['science','ieee','scatter','grid'])
plt.style.use(['science','grid'])

def norm(array):
    min_val = np.min(array,axis=1).reshape(-1,1)
    max_val = np.max(array,axis=1).reshape(-1,1)


    return array/max_val
   

def read_compute_norm(pathres,folder_root,sname):
    # Import results from the given path
    results = import_results(pathres)

    # initial_distribution, distribution, cliq, data, tf, dt, npart, records_p_dist

    # Determine the number of time steps


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
    # Compute mean particle concentration
    mean_p_c = (n_particle / vtot).reshape(-1,1)

    # Compute particle concentration variance
    par_var = np.sum(np.power(p_concentration - mean_p_c, 2) * full_volume, axis=1) / vtot
#    par_var = np.sqrt(par_var)

    # Compute mean concentration
    mean_c = np.sum(concentration_record * full_volume, axis=1) / vtot
    mean_c = mean_c.reshape(-1, 1)  # Reshape for broadcasting


    # Compute concentration variance
    var_c = np.sum(np.power(concentration_record - mean_c, 2) * full_volume, axis=1) / vtot

#    var_c = np.sqrt(var_c)
    # Normalize concentration record and particle concentration
    normalized_scalar_concentration = concentration_record / mean_c
    normalized_particle_concentration = p_concentration / mean_p_c

    # vtk_cma_mesh_path="/home/benjamin/Documenti/code/cpp/BIREM_new/out/sanofi/cma_mesh.vtu"
    # birem.vtk.mk_series(vtk_cma_mesh_path,"./results/",sname,results.t,
    # [full_volume,"liquid_volume"],[p_concentration,"particle_concentration"],[normalized_particle_concentration,"normalized_particle_concentration"],[normalized_scalar_concentration,"normalized_liquid_concentration"],[concentration_record,"liquid_concentration"])
    return normalized_particle_concentration,normalized_scalar_concentration,par_var/par_var[0],var_c/var_c[0],results.t

folder_root = "/home/benjamin/Documenti/code/cpp/biomc/cma_data/"

sname = ["data_50K","data_1M","data_5M","data_10M"]
# pathres = ['./results/50k.h5','./results/1M.h5','./results/5M.h5',"./results/10M.h5"]
#pathres = ['./results/50k.h5',"./results/1M.h5","./results/5M.h5"]
pathres = [ "./results/result_.h5"]

_,n_c ,_,pvc,t= read_compute_norm(pathres[0],folder_root,sname[0])
plt.semilogy(t,pvc,label="liquid")
for i in range(len(pathres)):
    n,n_c,par_var,_,t = read_compute_norm(pathres[i],folder_root,sname[i])
    plt.semilogy(t, par_var,label=sname[i])




plt.legend()
plt.title("Segregation index as a function of the time")
plt.ylabel(r"\[ \frac{\sigma(t)}{\sigma(t_{0})}\]")
plt.xlabel("time [s]")
plt.savefig("./results/mixing_variance2.svg",dpi=1500)
# plt.show()


