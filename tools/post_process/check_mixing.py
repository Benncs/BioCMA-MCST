from read_results import import_results
import numpy as np 
import matplotlib.pyplot as plt 

import birem.birem_generate
import birem 
import birem.vtk

def norm_distribution(raw_data):
  mind = np.min(raw_data,axis=0)
  maxd = np.max(raw_data,axis=0)
  return (raw_data-mind)/(maxd-mind)

def multipage(filename, figs=None, dpi=200):
  from matplotlib.backends.backend_pdf import PdfPages
  pp = PdfPages(filename)
  if figs is None:
    figs = [plt.figure(n) for n in plt.get_fignums()]
    for fig in figs:
      fig.savefig(pp, format='pdf')
      # plt.close(fig)
    pp.close()

def check_concentration_mixing(c_liq_records):
    cmax = np.max(c_liq_records[-1,:,0])
    normalized = data[-1,:,0]/c_liq_records[-1,0,0]
    diff = np.abs(normalized-1)*100
    err = np.where(diff>1)
    if len(err)>len(c_liq_records[0,:,0])//10:
      print("Non mixed yet")
      print(diff[err])
      print(err)

def plot_concentration_mixing(data):
    n_t = data.shape[0]
    t = np.linspace(0,tf,n_t)
    ind = [0,4,n_t//10,n_t//5,n_t-1]

    f = [ ]
    f.append(plt.figure())
    plt.title("Normalized Concentration")
    for i in ind:
        cmax = np.max(data[i,:,0])
        normalized = data[i,:,0]/cmax
        plt.plot(normalized,label=f"t={t[i]}s")
    plt.legend()
    plt.yscale('log')
    plt.ylabel("Normalized Concentration [-]")
    plt.xlabel("Compartment number")

 
    nc = len(distribution)
    ci = [0,nc//4,nc//2,nc-1]
    cref = data[-1,0,0]
    normalized = data[:,ci,0]/cref

    fig,axs = plt.subplots(1,2,sharey='row')
    f.append(fig)
    fig.suptitle("Normalized Concentration")
    axs[0].plot(data[0,:,0]/cref)
    axs[0].set_yscale('log')
    axs[1].set_yscale('log')
    axs[0].set_title("Initial condition") 
    axs[1].plot(data[-1,:,0]/cref)
    axs[1].set_title("Final time") 
    axs[1].set_xlabel("Compartment number")
    axs[0].set_xlabel("Compartment number")
    axs[0].set_ylabel("Normalized Concentration [-]")

    fig = plt.figure()
    f.append(fig)
    plt.title("Normalized Concentration")
    for i in ci:
        normalized_i = data[:,i,0]/cref
        plt.plot(t,normalized_i,"*-",label=f"{i}")
    plt.xlabel("Time [s]")
    plt.ylabel("Normalized Concentration [-]")
    plt.yscale('log')
    plt.legend()


    

    multipage("concentration_mixing.pdf",f,1500)

pathres = './results/result.h5'
initial_distribution,distribution,cliq,data,tf,dt,npart,records_d = import_results(pathres)



vtk_cma_mesh_path = "/home/benjamin/Documenti/cpp/BIREM_Project/out/sanofi/cma_mesh.vtu"
vtp_result = "./results/sanofi.vtu"

liquid_glucose_concentration_tmp = data[:,:,0]
liquid_glucose_concentration = np.zeros((liquid_glucose_concentration_tmp.shape[0]+1,liquid_glucose_concentration_tmp.shape[1]))
liquid_glucose_concentration[1:]=liquid_glucose_concentration_tmp
liquid_glucose_concentration[0]=liquid_glucose_concentration_tmp[0]

records_d_2 = np.zeros((records_d.shape[0]+1,records_d.shape[1]))
records_d_2[1:]=records_d[:,:,0]
records_d_2[0]=initial_distribution
records_d = records_d_2
# data = np.insert(data, 0, data[0])
# volumes = birem.vtk.read_scalar(vtk_cma_mesh_path,"liquid_volume")

liquid_volumes = birem.read_scalar("/home/benjamin/Documenti/cpp/BIREM_Project/out/sanofi/vofL.raw")

particle_concentration = records_d/liquid_volumes
n_t = particle_concentration.shape[0]

t = np.linspace(0,tf,n_t)


# c0_norm = data[0,:,0]
# cend_norm = data[-1,:,0]
# sc_c0 = birem.vtk.mk_scalar(c0_norm,"C0_norm")
# sc_cend = birem.vtk.mk_scalar(cend_norm,"Cend_norm")

# write json like this :
n = particle_concentration/np.max(particle_concentration,axis=1).reshape(-1,1) #norm_distribution(particle_concentration)
nc = liquid_glucose_concentration/np.max(liquid_glucose_concentration,axis=1).reshape(-1,1) #norm_distribution(liquid_glucose_concentration)


low = particle_concentration[:,particle_concentration[-1,:]<21600]
high = particle_concentration[:,particle_concentration[-1,:]>=21600]

print(np.where(particle_concentration[-1,:]>=21600))

nl  = low/np.max(low,axis=1).reshape(-1,1)

plt.plot(nl)
plt.show()

# birem.vtk.mk_series(vtk_cma_mesh_path,"./results/","sanofi",t,
# [particle_concentration,"particle_concentration"],
# [n,"normed_particle_concentration"],[liquid_glucose_concentration,"liquid_concentration"],
# [nc,"normed_liquid_concentration"])
    






# plt.plot(initial_distribution)

# plt.figure()
# plt.title("Normalized particule concentration")
# plt.plot(normalized_particle_concentration,label="Normalized particle concentration")
# plt.plot(initial_p_c,label="Initial normalized particle concentration")

# # plt.yscale("log")
# plt.legend()
# plt.show()

# VTK export 
# safe_points("6612_export",coordinates)
# birem.vtk.append_vtk(vtp_point_path,vtp_result,sc_c0,sc_cend,sc_p)

# plt.figure()
# plt.title("Normalized distribution (min-max normalization)")
# plt.plot(distribution,"*",label="particle distribution")
# plt.plot(volumes,label="volume distribution")
# plt.legend()
# plt.show()