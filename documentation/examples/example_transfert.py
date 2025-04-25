"""
Minimal working program to use BioMC API V0.5.

This program demonstrates a basic example of how to perform a simulation using BioMC from a Python program.
It focuses on:
- Creating a handle for shared execution (No MPI)
- Performing basic simulation settings and running it
- Loading and running a basic simulation

Note that the linked library is the shared one, this will not work linking the distributed one.

Usage:
    Run this script to perform a simulation using BioMC API V0.5.
"""

# Import the handle module
import handle_module
import numpy as np 
import os 
import cmtool
import tempfile
import biomc_pp
import matplotlib.pyplot as plt 
import shutil
import plotly.graph_objs as go
def mk_flowmap(dest,name):
    eps_turb = 0.01
    outfolder = f"{dest}/{name}"
    os.makedirs(outfolder,exist_ok=True)
    _case = cmtool.generate.generate_0D_set_from_fraction(outfolder,20e-3,5/100,eps_turb)
    cmtool.cma_case.mk_write_cma_case(f"{outfolder}/cma_case",[1,0,0],description=f"{name}",cma_path=_case)
    return outfolder

def run(params,out,name,cma_path): 
    handle = handle_module.init_simulation(out,name,cma_path,params)
    liquid_concentration_0 = np.zeros((1,2))  
    gas_concentration_0 = np.zeros((1,2)) 
    gas_concentration_0[:,1]=300e-3
    handle_module.set_initial_concentrations(handle,liquid_concentration_0,gas_concentration_0)
    handle_module.register_model_name(handle, "None")
    rc = handle_module.apply(handle, False)
    if not rc[0]:
        print(rc[1])
        return -1
    rc = handle_module.exec(handle)

def show_result(dirname):
    post_process = biomc_pp.PostProcess("example_transfert",dirname)
    time = biomc_pp.check_time_unit(post_process)
    C_o_l = post_process.get_spatial_average_concentration(1,biomc_pp.Phase.Liquid)
    C_o_g = post_process.get_spatial_average_concentration(1,biomc_pp.Phase.Gas)
    fig = go.Figure()

    fig.add_trace(go.Scatter(
        x=time,
        y=C_o_l,
        mode='lines',
        name='Liquid Concentration',
        line=dict(color='blue'),
        yaxis='y1'
    ))

    fig.add_trace(go.Scatter(
        x=time,
        y=C_o_g,
        mode='lines',
        name='Gas Concentration',
        line=dict(color='red'),
        yaxis='y2'
    ))

    fig.update_layout(
        title='Spatial Average Concentration',
        xaxis=dict(title='Time'),
        yaxis=dict(
            title='Liquid Concentration',
        ),
        yaxis2=dict(
            title='Gas Concentration',
            overlaying='y',
            side='right',
        ),
        legend=dict(x=0.5, y=1.1, orientation='h', xanchor='center'),
    )

    fig.show()

if __name__ == "__main__":

    params = handle_module.make_params(
        biomass_initial_concentration=0.1,
        final_time=100,
        delta_time=0.01,
        number_particle=1,
        number_exported_result=50,
    )

    with tempfile.TemporaryDirectory() as tmpdirname:
        cma_path = mk_flowmap(tmpdirname,"0d_gas")
        run(params,tmpdirname,"example_transfert",cma_path+"/")
        show_result(tmpdirname)
        shutil.rmtree(tmpdirname)


        
