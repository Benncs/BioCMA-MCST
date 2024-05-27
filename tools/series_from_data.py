import birem 
import numpy as np 
import birem.vtk 
from read_results import * 
import matplotlib.pyplot as plt 


def norm(array):
    min_val = np.min(array,axis=1).reshape(-1,1)
    max_val = np.max(array,axis=1).reshape(-1,1)

    # S'assurer que la différence n'est pas zéro pour éviter la division par zéro
    # if max_val - min_val != 0:
    # return (array - min_val) / (max_val - min_val)
    return array/max_val
    # return array/max_val
    # else:
    #     return array - min_val
      




def read_compute_norm(pathres,folder_root,sname):
    initial_distribution,distribution,cliq,data,tf,dt,npart,records_d = import_results(pathres)

    n_t = records_d.shape[0]
    read_volume = np.zeros((14, records_d.shape[1]))
    t = np.linspace(0, tf, n_t)

    t_per_flowmap = 0.0286
    n_per_flowmap = int(t_per_flowmap / dt)

    # Read volume data from files and store in read_volume
    for i in range(14):
        folder = f"/{folder_root}/bench/i_{i+1}"
        raw_volume = birem.read_scalar(folder + "/vofL.raw")
        read_volume[i] = raw_volume

    records_d = records_d.reshape(records_d.shape[0], records_d.shape[1])

    p_concentration = np.zeros((records_d.shape[0],records_d.shape[1]))
    full_volume = np.zeros_like(p_concentration)

    full_volume[0] = np.copy(read_volume[0])
    p_concentration[0] = records_d[0] / read_volume[0]

    concentration_record = data[:,:,0]

    for i in range(1, n_t):
        index = (i // n_per_flowmap) % 14
        full_volume[i] = np.copy(read_volume[index])
        p_concentration[i] = records_d[i] / read_volume[index]

    print(np.sum(records_d[0]))
    n = norm(p_concentration)# / np.max(p_concentration, axis=1).reshape(-1, 1)
    n_c = norm(concentration_record) #/ np.max(concentration_record, axis=1).reshape(-1, 1)


    print(np.std(p_concentration[-1])/np.mean(p_concentration[-1])*100)
    print(np.std(concentration_record[-1])/np.mean(concentration_record[-1])*100)
    vtk_cma_mesh_path="/home/benjamin/Documenti/code/cpp/BIREM_new/out/cma_mesh.vtu"
    birem.vtk.mk_series(vtk_cma_mesh_path,"./results/",sname,t,
    [full_volume,"liquid_volume"],[p_concentration,"particle_concentration"],[n,"normalized_particle_concentration"],[n_c,"normalized_liquid_concentration"],[concentration_record,"liquid_concentration"])
    return n,n_c

# s_core = (p_concentration[-1]-np.mean(p_concentration[-1]))/(np.std(p_concentration[-1]))

#pathres = './results/50k.h5'
folder_root = "/home/benjamin/Documenti/code/cpp/biomc/cma_data/"

sname = ["data_50K","data_1M","data_5M","data_10M"]
pathres = ['./results/50k.h5','./results/1M.h5','./results/5M.h5',"./results/10M.h5"]

for i in range(len(pathres)):
    n,n_c = read_compute_norm(pathres[i],folder_root,sname[i])

    plt.plot(n[-1],label=pathres[i])
    print("max",np.max(n,axis=1))
# plt.plot(n[-1],color='black')
# plt.yscale('log')
_,n_c = read_compute_norm(pathres[0],folder_root,sname[0])
plt.plot(n_c[-1],label="liquid")
plt.legend()
plt.show()


