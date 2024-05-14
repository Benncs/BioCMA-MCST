import numpy as np 
import h5py 

def import_results(path):
  with h5py.File(path, 'r') as file:
      # Create a group for attributes
      distribution = np.array(file["final_results/distribution"])
      cliq = np.array(file["final_results/concentrations/liquid"])
      data = np.array(file["final_results/concentrations/records"])
      tf = np.array(file['initial_parameters/final_time'])
      dt = np.array(file['initial_parameters/delta_time'])
  return distribution,cliq,data,tf,dt