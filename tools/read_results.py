import numpy as np 
import h5py 

def import_results(path):
  with h5py.File(path, 'r') as file:
      # Create a group for attributes
      distribution = np.array(file.get("final_results/distribution"))
      initial_distribution = np.array(file.get("initial_parameters/particle_distribution"))
      cliq = np.array(file.get("final_results/concentrations/liquid"))
      data = np.array(file.get("final_results/concentrations/records"))
      tf = np.array(file.get('initial_parameters/final_time'))
      dt = np.array(file.get('initial_parameters/delta_time'))
      npart  = np.array(file['initial_parameters/number_particles_0'])
      
      records_distribution = np.array(file.get("final_results/concentrations/records_distribution"))
  return initial_distribution,distribution,cliq,data,tf,dt,npart,records_distribution