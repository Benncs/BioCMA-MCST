
import numpy as np
# from biomc import make_initial_concentration

import h5py 
import numpy as np
from typing import Optional


def make_initial_concentration(dest:str,liq_concentration:np.ndarray,gas_concentration:Optional[np.ndarray]=None)->None:
    data_liq = np.asfortranarray(liq_concentration)
    data_gas = np.zeros(())
    if gas_concentration is not None:
        data_gas = np.asfortranarray(gas_concentration)

    with h5py.File(dest, "w") as file: 
        file.create_dataset("initial_liquid", data=data_liq)
        if(gas_concentration is not None):
            file.create_dataset("initial_gas", data=data_gas)



def init_0d_1s():
    liquid_0d = np.zeros((1,2))
    liquid_0d[0,0]=0
    make_initial_concentration("./cma_data/0d_init.h5",liquid_0d)

def init_0d_4s():
    liquid_0d = np.zeros((1,4))
    gas_0d = np.zeros((1,4))
    liquid_0d[0,0]=5
    liquid_0d[0,1]=0.

    gas_0d[0,0]=0.
    gas_0d[0,1]=300e-3

    make_initial_concentration("./cma_data/0d_4s_init.h5",liquid_0d,gas_0d)

def init_n14_1s():
    liquid_n14 = np.zeros((500,1))
    liquid_n14[:,0]=0.
    liquid_n14[0,0]=5.
    make_initial_concentration("./cma_data/n14_init.h5",liquid_n14)

def init_n14_4s():
    liquid_sanofi = np.zeros((500,4))
    gas_sanofi = np.zeros((500,4))

    gas_sanofi[:,0]=0 #glucose g/l
    gas_sanofi[:,1]=0.21 #o2 g/l

    liquid_sanofi[:,0]=5 #glucose g/l
    liquid_sanofi[:,1]=0 #o2 g/l

    make_initial_concentration("./cma_data/n14_4s_init.h5",liquid_sanofi,gas_sanofi)

def sanofi_1s():
    liquid_sanofi = np.zeros((432,4))
    gas_sanofi = np.zeros((432,4))


    gas_sanofi[:,0]=0 #glucose g/l
    gas_sanofi[:,1]=0.21 #o2 g/l

    liquid_sanofi[:,0]=5 #glucose g/l
    liquid_sanofi[:,1]=1e-3 #o2 g/l


    make_initial_concentration("./cma_data/sanofi_init.h5",liquid_sanofi,gas_sanofi)

if __name__=="__main__":
    
    # liquid_0d = np.zeros((1,2))
    # liquid_n50 = np.zeros((50,1))

    # liquid_n50[:,0]=0.
    
    # make_initial_concentration("./cma_data/n50_init.h5",liquid_n50)
    
    init_0d_4s()
    init_0d_1s()
    init_n14_1s()
    sanofi_1s()
    init_n14_4s()