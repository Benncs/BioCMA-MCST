import matplotlib.pyplot as plt
import numpy as np 
import scienceplots
plt.style.use(['science','ieee','no-latex'])
N0 = np.array([200000, 200000, 200000, 200000, 200000, 10000, 10000, 10000, 10000])
equivalent_mu = np.array([0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.77, 0.8])
Q = np.array([5.56E-07, 1.11E-06, 1.67E-06, 2.22E-06, 2.78E-06, 3.33E-06, 3.89E-06, 4.27777777777778E-06, 42.7777777777778])
Tau = np.array([3.60E+04, 1.80E+04, 1.20E+04, 9.00E+03, 7.20E+03, 6.00E+03, 5.14E+03, 4.68E+03, 4.68E-04])
D = np.array([2.78E-05, 5.56E-05, 8.33E-05, 1.11E-04, 1.39E-04, 1.67E-04, 1.94E-04, 2.14E-04, 2.22E-04])
N = np.array([None, 3.77E+05, 3.56E+05, 347926, 3.42E+05, 1.72E+04, 16217, 3987, 1219])
X = np.array([None, 2.4, 2.46247611104735, 2.47962229210562, 2.45, 2.45340750554035, 2.46247611104735, 0.6, 0.176031616998231])
mu = np.array([None, 0.00005897, 0.000085409, 0.00011229, 0.00014023, 0.00016862, 0.00019588, 0.00021333, 0.00021333])
S = np.array([None, 0.0036916, 0.0065321, 0.010949, 0.018827, 0.037038, 0.11538, 3.8309, 4.6451])



D_valid = D[1:]*3600
X_valid = np.array(X[1:]) 
S_valid = np.array(S[1:])  

plt.figure(figsize=(10, 6))
plt.plot(D_valid, X_valid, '-*', label='X')
plt.plot(D_valid, S_valid, '-*', label='S')

plt.xlabel(r'D [$h^{1}$]')
plt.ylabel('Concentration [G/L]')
plt.title('Biomass and Glucose concentration at steady state according to dillution rate')
plt.grid(True)
plt.legend()
plt.savefig("tmp")


