import birem 
import numpy as np 
import birem.vtk 
from read_results import * 
import matplotlib.pyplot as plt 


def norm(array):
    min_val = np.min(array)
    max_val = np.max(array)

    # S'assurer que la différence n'est pas zéro pour éviter la division par zéro
    if max_val - min_val != 0:
        return (array - min_val) / (max_val - min_val)
    else:
        return array - min_val
      


pathres = './results/result.h5'
folder_root = "/home/benjamin/Documenti/cpp/BioCMA-MCST/cma_data/"
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



for i in range(1, n_t):
    index = (i // n_per_flowmap) % 14
    full_volume[i] = np.copy(read_volume[index])
    p_concentration[i] = records_d[i] / read_volume[index]

print(np.sum(records_d[0]))
n = p_concentration / np.max(p_concentration, axis=1).reshape(-1, 1)

# plt.plot(norm(p_concentration[-1]))
# # plt.plot(n[-1],color='black')
# # plt.yscale('log')
# plt.show()

vtk_cma_mesh_path="/home/benjamin/Documenti/cpp/BIREM_Project/out/data/cma_mesh.vtu"
birem.vtk.mk_series(vtk_cma_mesh_path,"./results/","data_5m",t,
[full_volume,"liquid_volume"],[p_concentration,"particle_concentration"],[n,"normed_particle_concentration"])
