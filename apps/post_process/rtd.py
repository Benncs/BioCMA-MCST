import numpy as np

from . import FIGURE_TYPE, get_time
from .read_results import Results
import matplotlib.pyplot as plt 

def str_rtd(t: np.ndarray, tau):
    return 1 / tau * np.exp(-t / tau) # This work only for cstr

def get_tau_from_rtd(rtd):
    return 1/rtd[0] # This work only for cstr


def get_particle_rtd(given_flow: float, results):
    num_bins = 100
    probes = np.zeros(1,)

    for i in results.partial:
        if i.probes is not None:  # Ensure probes are not None
            probes = np.concatenate((probes, i.probes))


    c, e = np.histogram(probes, bins=num_bins,density=True)
    return c ,e




def get_scalar_rtd(dest_root:str,given_flow: float, results,is_str:bool=False) -> None:
    
    
    # Get the concentration everywhere at t=0
    c0 = results.main.concentrations_liquid[0, :, 0]
    # Get the concentration everywhere for everytimestep but we substrat initial concentration to only have Delta C
    delta_c = results.main.concentrations_liquid[:, :, 0] - c0


    step_concentration = 5

    # Cumulative probability 
    rtd_f = delta_c[:, -1] / step_concentration  
    # By defintion F(t)=integral 0 to t (E(t))
    e_rtd = np.diff(rtd_f) / np.diff(results.time)
    
    

    counts_normalized,bin_edges= get_particle_rtd(given_flow,results)
    V = np.sum(results.main.volume_liquid, axis=1)[0]
    Q = given_flow
    tau_analytical = V / Q
    print(f"Analytical tau :{tau_analytical}s")
    
    if(is_str):
        
        n_t_tau = np.where(np.abs((results.time - tau_analytical)) < 100)
        print(np.sum(results.total_repartion,axis=1)[n_t_tau[0]].shape)
        reduction = (
            100 - np.sum(results.total_repartion,axis=1)[n_t_tau] / np.sum(results.total_repartion,axis=1)[0] * 100
        )
        rtd_str_analytical = str_rtd(results.time, tau_analytical)
        
        simulated_tau = get_tau_from_rtd(e_rtd) 
        
        print(f"Tau from scalar simulation {simulated_tau/3600}h")
        print(f"Particule number reduction after 1 tau (expected 63%): {reduction}")


    plt.figure()
    plt.plot(results.time,rtd_f,'black',label='simulated cumulative rtd')
    plt.legend()
    plt.xlabel(f"Time [{get_time()}]")
    plt.ylabel("F(t)")
    plt.savefig(f"{dest_root}/rtd_cumsum{FIGURE_TYPE}")

    plt.figure()
    #start from 1 because of np.diff
    plt.plot(results.time[1:],e_rtd,'--',color='red',label='simulation')
    if(is_str):
        plt.plot(results.time,rtd_str_analytical,color="black",label="theoretical rtd")

    plt.xlabel(f"Time [{get_time()}]")
    plt.ylabel("E(t)")
    plt.legend()
    plt.savefig(f"{dest_root}/rtd{FIGURE_TYPE}")
  

    plt.figure()
    plt.bar(
        bin_edges[:-1],
        counts_normalized,
        width=np.diff(bin_edges),
        edgecolor="black",
        alpha=0.7,
        color="blue",
        label="Particles"
    )
    plt.xlabel(f"Time [{get_time()}]")
    plt.plot(results.time[1:],e_rtd,'--',color='red',label='Scalar')
    plt.ylabel("E(t)")
    plt.legend()
    plt.savefig(f"{dest_root}/rtd_normalized{FIGURE_TYPE}")
