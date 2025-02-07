import biomc_pp
import matplotlib.pyplot as plt 
import numpy as np 
pp = biomc_pp.PostProcess('test_sanofi',"/home-local/casale/Documents/thesis/code/BioCMA-MCST/results/")



nu1= pp.get_population_mean("nu1",pp.max_n_export_bio-1)
nu2= pp.get_population_mean("nu1",pp.max_n_export_bio-1)
nu_1_t = pp.get_time_population_mean("nu1")
nu_2_t = pp.get_time_population_mean("nu2")
m_t= pp.get_time_population_mean("mass")
m= pp.get_population_mean("mass",pp.max_n_export_bio-1)
mu = (nu1+nu2)/m




v=299#20e-3

# v=20e-3

time = biomc_pp.check_time_unit(pp)

# X = pp.get_biomass_concentration()
weigth=529110468072872.7
X = np.array(
        [
            np.sum(weigth * np.array(pp.get_properties("mass", i)))
            for i in range(pp.n_export)
        ]
)/ v



# Getting the data
S = pp.get_spatial_average_concentration(0, biomc_pp.Phase.Liquid)
O2 = pp.get_spatial_average_concentration(1, biomc_pp.Phase.Liquid)
A = pp.get_spatial_average_concentration(2, biomc_pp.Phase.Liquid)

fig, (ax1, ax2, ax3,ax4) = plt.subplots(4, 1, figsize=(10, 15))

ax1.plot(time, X, color="red", label="X")
ax1.plot(time, S, color="black", label="S")
ax1.set_xlabel('Time')
ax1.set_ylabel('X and S')
ax1.legend(loc='upper left')

ax11 = ax1.twinx()  
ax11.plot(time, O2, "--", label="O2", color="blue")
ax11.set_ylabel('O2')
ax11.legend(loc='upper right')

ax1.grid(True)

O2_gas = pp.get_spatial_average_concentration(1, biomc_pp.Phase.Gas)
ax2.plot(time, O2_gas, label="O2 Gas")
ax2.set_xlabel('Time')
ax2.set_ylabel('O2 Gas Concentration')
ax2.legend(loc="upper left")
ax2.grid(True)

ax3.plot(time, A, label="Acetate")
ax3.set_xlabel('Time')
ax3.set_ylabel('Acetate Concentration')
ax3.legend(loc="upper left")
ax3.grid(True)

ax4.plot(time, pp.get_growth_in_number(), label="N")
ax4.set_xlabel('Time')
ax4.set_xlabel('Growth in Number')
ax4.legend()
ax4.grid(True)
plt.tight_layout()



m0 = X[0] * v
mf = X[-1] * v
y_xs_app = (X[-1] - X[0]) / (S[0] - S[-1])
print(f"y_xs_app: {y_xs_app}")

dmdt = (mf - m0) / time[-1] * 3600

fig, ax1 = plt.subplots()
ax1.semilogy(time, nu_1_t, color="black", label="nu1")
ax1.set_xlabel('Time')
ax1.set_ylabel('nu1')

ax2 = ax1.twinx()
ax2.semilogy(time, nu_2_t, label="nu2")
ax2.set_ylabel('nu2')

fig.legend(loc='upper right')
ax1.grid(True)

plt.figure()
plt.plot(time[1:], pp.get_time_population_mean("mu")[1:] * 3600, color="black", label='mu')
plt.plot(time, (nu_1_t + nu_2_t) / m_t * 3600, label="mu_p_avg")
plt.xlabel('Time')
plt.ylabel('mu and mu_p_avg')
plt.legend()
plt.grid(True)

# Set up the seventh plot (a_pts and a_permease vs. Time)
plt.figure()
plt.plot(time, pp.get_time_population_mean("a_pts"), "-*", color="black", label='a_pts')
plt.plot(time, pp.get_time_population_mean("a_permease"), label="a_permease")
plt.xlabel('Time')
plt.ylabel('a_pts and a_permease')
plt.legend()
plt.grid(True)

# Show all the plots
plt.show()

