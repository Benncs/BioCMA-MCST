import numpy as np
from .read_results import Results
import matplotlib.pyplot as plt 

def str_rtd(t: np.ndarray, tau):
    return 1 / tau * np.exp(-t / tau) # This work only for cstr

def get_tau_from_rtd(rtd):
    return 1/rtd[0] # This work only for cstr


def get_particle_rtd(given_flow: float, results):
    pass


def get_scalar_rtd(given_flow: float, results) -> None:
    V = np.sum(results.main.volume_liquid, axis=1)[0]
    Q = given_flow
    tau_analytical = V / Q

    # Get the concentration everywhere at t=0
    c0 = results.main.concentration_liquid[0, :, 0]
    # Get the concentration everywhere for everytimestep but we substrat initial concentration to only have Delta C
    delta_c = results.main.concentration_liquid[:, :, 0] - c0


    step_concentration = 5

    # Cumulative probability 
    rtd_f = delta_c[:, -1] / step_concentration  

    n_t_tau = np.where(np.abs((results.time - tau_analytical)) < 100)
    reduction = (
        100 - results.total_repartion[n_t_tau] / results.total_repartion[0] * 100
    )

    rtd_str_analytical = str_rtd(results.time, tau_analytical)

    # By defintion F(t)=integral 0 to t (E(t))
    e_rtd = np.diff(rtd_f) / np.diff(results.time)
    simulated_tau = get_tau_from_rtd(e_rtd) 

    print(f"Analytical tau :{tau_analytical/3600}h")
    print(f"Tau from scalar simulation {simulated_tau/3600}h")

    #start from 1 because of np.diff
    plt.plot(results.time[1:],e_rtd,'--',color='red',label='simulation')
    # plt.plot(t,rtd_str,'black',label='rtd')
    plt.legend()
    plt.savefig("./tmp")
