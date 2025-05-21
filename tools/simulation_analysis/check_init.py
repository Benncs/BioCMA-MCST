import biomc_pp
import matplotlib.pyplot as plt
import numpy as np


pp = biomc_pp.PostProcess("debug","./results/")

l0 = pp.get_properties("length",0)

print(pp.estimate(biomc_pp.Estimator.MonteCarlo,"length",0)*1e6)

print(np.std(l0))

plt.hist(l0,bins=500)

plt.show()