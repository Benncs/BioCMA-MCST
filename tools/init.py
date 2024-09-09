
import numpy as np
from biomc import make_initial_concentration

if __name__=="__main__":
    # c =np.random.random((500,2))
    liquid_0d = np.zeros((1,2))
    liquid_n14 = np.zeros((500,2))
    liquid_0d[0,0]=1.
    liquid_n14[:,0]=1./500.

    make_initial_concentration("./cma_data/0d_init.h5",liquid_0d)
    make_initial_concentration("./cma_data/n14_init.h5",liquid_n14)