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
