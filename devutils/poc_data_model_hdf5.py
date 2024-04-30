import h5py
import numpy as np

# Create a new HDF5 file
with h5py.File('output.h5', 'w') as file:
    # Create a group for attributes
    
    # Define attributes
    file.attrs['creation_date'] = '2024-05-01'
    file.attrs['author'] = 'John Doe'
    file.attrs['description'] = 'Simulation output'
    
    # Create a group for data
    data_group = file.create_group('data')
    
    # Define data fields
    data_group.attrs['n_particles'] = np.uint64(1000)
    data_group.attrs['n_compartments'] = np.uint64(5)
    data_group.attrs['final_time'] = np.float64(10.0)
    data_group.attrs['delta_t'] = np.float64(0.1)
    
    # Create datasets for concentration matrices
    liquid_concentration = np.random.rand(3, 5)  # Example liquid concentration matrix
    gas_concentration = np.random.rand(3, 5)     # Example gas concentration matrix
    data_group.create_dataset('concentrations/liquid', data=liquid_concentration)
    data_group.create_dataset('concentrations/gas', data=gas_concentration)
    
    # Create datasets for integrated particle states
    userdefined1 = np.random.rand(5)  # Example data for userdefined1
    userdefined2 = np.random.rand(5)  # Example data for userdefined2
    data_group.create_dataset('integrated_particle_states/userdefined1', data=userdefined1)
    data_group.create_dataset('integrated_particle_states/userdefined2', data=userdefined2)
