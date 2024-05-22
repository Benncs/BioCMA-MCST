from read_results import import_results
import numpy as np 
import matplotlib.pyplot as plt 

vliq_bin = [4.44991e-05,4.44771e-05,4.29673e-05,3.33303e-05,4.69845e-05,5.24179e-05,5.24129e-05,5.16096e-05,5.13444e-05,5.1332e-05,4.40372e-05,4.39985e-05,4.19968e-05,3.16794e-05,4.46853e-05,4.95365e-05,4.95366e-05,4.95319e-05,4.95274e-05,4.95595e-05,4.02903e-05,4.02376e-05,4.07098e-05,3.26117e-05,4.51157e-05,4.98775e-05,4.98801e-05,4.98765e-05,4.98724e-05,4.98739e-05,4.72627e-05,4.72244e-05,4.42903e-05,3.2928e-05,4.75978e-05,5.33965e-05,5.33924e-05,5.33917e-05,5.33957e-05,5.3388e-05,4.03472e-05,4.03112e-05,4.05924e-05,3.20896e-05,4.83684e-05,5.47482e-05,5.47393e-05,5.50453e-05,5.50945e-05,5.50774e-05,4.29949e-05,4.29401e-05,4.22171e-05,3.32506e-05,4.64913e-05,5.1891e-05,5.18898e-05,5.13981e-05,5.12372e-05,5.12468e-05,4.26268e-05,4.2593e-05,4.2045e-05,3.29148e-05,4.41489e-05,4.75716e-05,4.75766e-05,4.89156e-05,4.9104e-05,4.91111e-05,4.28321e-05,4.28009e-05,4.16638e-05,3.24315e-05,4.94113e-05,5.68563e-05,5.68622e-05,5.54631e-05,5.49903e-05,5.4984e-05,4.52837e-05,4.52403e-05,4.32924e-05,3.2751e-05,4.53496e-05,4.97791e-05,4.97821e-05,4.97899e-05,4.97909e-05,4.97767e-05,4.58668e-05,4.58455e-05,4.34052e-05,3.24781e-05,4.52243e-05,4.92544e-05,4.92557e-05,5.14376e-05,5.17547e-05,5.17151e-05,0.000109936,0.000109817,0.000116816,0.000102211,0.000118929,0.000121736,0.000121783,0.000121796,0.000121789,0.000121769,0.000107905,0.000107827,0.000117272,0.000105409,0.000120246,0.000122991,0.000122991,0.000122987,0.000122983,0.000122946,0.00010887,0.000108756,0.000115925,0.000101011,0.000120444,0.000124564,0.000124559,0.000124556,0.000124554,0.000124506,0.000110155,0.00011005,0.000115886,0.00010186,0.000119535,0.000124357,0.00012435,0.000124353,0.000124353,0.000124252,0.000107881,0.000107756,0.000116656,0.000102244,0.000121199,0.000124617,0.00012463,0.000124631,0.000124622,0.000124649,0.000108535,0.000108433,0.000116202,0.000102086,0.0001191,0.000122148,0.000122198,0.000122202,0.000122207,0.000122172,0.000110767,0.000110674,0.000118331,0.00010509,0.000122014,0.000125699,0.000125761,0.000125761,0.000125762,0.00012571,0.000111557,0.000111447,0.000117027,0.000101039,0.00011874,0.000121939,0.000121947,0.000122509,0.000122513,0.000121919,0.00010775,0.00010765,0.000114721,0.000101695,0.000120779,0.000126265,0.00012634,0.000126192,0.000126195,0.000126318,0.000107952,0.000107853,0.000116956,0.000102673,0.000118764,0.000120518,0.000120498,0.00012055,0.000120551,0.000120524,0.000165772,0.000163169,0.000175913,0.000160075,0.000176351,0.000180139,0.000179406,0.000178233,0.000178197,0.000178161,0.000167263,0.000167125,0.00017842,0.000159094,0.000175341,0.000178255,0.00017826,0.000178246,0.000178239,0.000178218,0.000161669,0.000161528,0.000175083,0.000158991,0.000174486,0.000177913,0.000177903,0.000177896,0.000177908,0.000177893,0.000165369,0.000165221,0.000176828,0.00015847,0.000174092,0.00017715,0.000177135,0.00017715,0.000177158,0.000177115,0.000162007,0.000164121,0.000177508,0.000159136,0.000174693,0.000177702,0.000177711,0.000178222,0.000179663,0.000179637,0.000168819,0.000166768,0.000178076,0.000159998,0.00017358,0.000175161,0.00017516,0.00017464,0.000172534,0.00017251,0.000165639,0.000165478,0.000177834,0.000160482,0.00017672,0.000179397,0.000179423,0.000180903,0.000180562,0.000180542,0.000163664,0.000163518,0.000176092,0.000159024,0.000171218,0.000173343,0.000173335,0.000172666,0.000172021,0.000171998,0.000166294,0.00016614,0.000177123,0.000158297,0.000175292,0.000179138,0.000179126,0.000179155,0.000179164,0.000179132,0.000165356,0.000166931,0.000178664,0.000158815,0.000172258,0.000173734,0.000175042,0.000176332,0.00017635,0.000176338,0.000270449,0.000267724,0.000296062,0.000279779,0.000282042,0.000288022,0.000286829,0.00028219,0.000281991,0.000281203,0.000276901,0.000276795,0.000305855,0.000287329,0.000288806,0.000294908,0.000294885,0.000294876,0.000294869,0.000294846,0.000257614,0.000257529,0.000284526,0.000267278,0.000268628,0.000274304,0.000274301,0.000274299,0.000274307,0.000274288,0.000276931,0.000276837,0.000305876,0.000287309,0.000288765,0.00029486,0.000294856,0.000294862,0.000294865,0.000294839,0.000265021,0.000267737,0.000295595,0.000275066,0.000275707,0.000280819,0.000280312,0.000282414,0.00028797,0.000287931,0.000269902,0.000267556,0.000295968,0.000279882,0.000282041,0.000288008,0.000288662,0.000287932,0.000282152,0.000282539,0.000276937,0.000276834,0.000305868,0.000287314,0.000288766,0.000294874,0.000294903,0.000294886,0.000294862,0.000294844,0.000257624,0.000257534,0.000284547,0.000267282,0.000268642,0.000274321,0.000274305,0.000274302,0.00027431,0.000274289,0.000276922,0.000276817,0.000305864,0.00028731,0.000288771,0.000294866,0.000294857,0.000294863,0.000294866,0.00029484,0.000264063,0.000267609,0.000295478,0.000275291,0.000275711,0.000280824,0.000283234,0.00028799,0.000288326,0.000289266,0.00025673,0.000256729,0.000284193,0.000267282,0.000267634,0.00027252,0.000272527,0.000272541,0.000272543,0.000272625,0.000286973,0.000286392,0.000318112,0.000298165,0.000299746,0.000304016,0.000303997,0.000304016,0.000303975,0.00030478,0.000251621,0.000251619,0.000278467,0.000261897,0.000262289,0.000267112,0.000267117,0.000267129,0.000267112,0.000267214,0.000286979,0.000286398,0.000318136,0.000298189,0.000299741,0.000303971,0.000303981,0.000304016,0.000303976,0.00030478,0.000256737,0.000256736,0.000284192,0.000267276,0.000267633,0.000272519,0.000272489,0.000272495,0.000272509,0.000272622,0.000256731,0.000256734,0.000284195,0.000267282,0.000267638,0.000272525,0.000272496,0.000272501,0.000272522,0.000272634,0.000286979,0.000286393,0.000318119,0.000298189,0.000299739,0.000303971,0.000304003,0.000304061,0.00030401,0.000304782,0.000251617,0.000251616,0.000278468,0.000261851,0.000262283,0.000267134,0.000267137,0.000267122,0.000267105,0.000267207,0.000286972,0.000286393,0.000318123,0.000298216,0.000299749,0.000303971,0.000303981,0.000304016,0.000303976,0.00030478,0.000256737,0.000256736,0.000284192,0.000267276,0.000267633,0.000272519,0.000272489,0.000272495,0.000272509,0.00027264]

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
# check_concentration_mixing(data)
# plot_concentration_mixing(data)

# particle_distribution = norm_distribution(distribution)
# volume_distribution = norm_distribution(volumes)

# plt.show()
# exit(0)

# def write_json_series(destination,name,file_entries):
#   import json
#   json_content = {
#     "file-series-version": "1.0",
#     "files": file_entries
# }

#   # Write the JSON file
#   json_file_path = f"{destination}/{name}.vtu.series"
#   with open(json_file_path, 'w') as json_file:
#       json.dump(json_content, json_file, indent=4)

# def mk_series(filepath,destination,series_name,t,*args):
#   file_entries = []
#   n_t = len(t)
#   for i in range(n_t):
#     vtp_result = f"./{destination}/{series_name}_{i}.vtu"
#     scalar =[]
#     for s in args:
#       data = s[0]
#       name = s[1]
#       scalar.append(birem.vtk.mk_scalar(data[i],name))
#     birem.vtk.append_scalar(filepath, vtp_result,*scalar)
#     file_entries.append({"name": f"{series_name}_{i}.vtu", "time": t[i]})
#   write_json_series(destination,series_name,file_entries)




import birem.birem_generate
import birem 
import birem.vtk

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