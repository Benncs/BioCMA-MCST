
import numpy as np
from biomc import make_initial_concentration

if __name__=="__main__":
    
    liquid_0d = np.zeros((1,2))
    liquid_n14 = np.zeros((500,1))
    liquid_n50 = np.zeros((50,1))
    liquid_0d[0,0]=1
    liquid_n14[:,0]=0.
    liquid_n50[:,0]=0.
    # make_initial_concentration("./cma_data/0d_init.h5",liquid_0d)
    
    make_initial_concentration("./cma_data/n14_init.h5",liquid_n14)
    make_initial_concentration("./cma_data/n50_init.h5",liquid_n50)