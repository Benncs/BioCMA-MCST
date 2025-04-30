import numpy as np 
import matplotlib.pyplot as plt 
import biomc_pp

# def get_particle_rtd(given_flow: float, results):
#     num_bins = 100
#     probes = np.zeros(1,)

#     for i in results.partial:
#         if i.probes is not None:  # Ensure probes are not None
#             probes = np.concatenate((probes, i.probes))


#     c, e = np.histogram(probes, bins=num_bins,density=True)
#     return c ,e


pp = biomc_pp.PostProcess("rtd_sanofi","/home-local/casale/Documents/code/poc/results/")
time = biomc_pp.check_time_unit(pp)
probes = pp.get_probes()/3600
print(probes)

tallies = pp.get_tallies()

exit_tallies = tallies[:,3]

# plt.hist(exit_tallies,bins=100,density=True)
# plt.plot(time,exit_tallies)
# plt.bar(exit_tallies)
plt.hist(probes,density=True,bins=50,    edgecolor="black",
        alpha=0.7,
        color="blue",
        label="Particles")



plt.show()


# import h5py
# import numpy as np
# import matplotlib.pyplot as plt

# # Open the HDF5 file
# file_path = '/home-local/casale/Documents/code/poc/results/rtd_sanofi/rtd_sanofi_partial_0.h5'  # Replace with your HDF5 file path
# with h5py.File(file_path, 'r') as file:
#     # Read the dataset called 'probes' into a NumPy array
#     probes_dataset = file['probes'][:]

# # Create a histogram of the dataset
# plt.hist(probes_dataset, bins=100, edgecolor='black')  # Adjust the number of bins as needed
# plt.title('Histogram of Probes Dataset')
# plt.xlabel('Value')
# plt.ylabel('Frequency')
# plt.show()
