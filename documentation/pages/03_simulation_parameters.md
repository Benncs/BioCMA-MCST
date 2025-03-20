\page simparam  Simulation Configuration

## Scalar initialisation

There are different ways to initialize the initial concentration, but one simple method is to generate a .h5 file using the tools/init.py script.

A template fonction to add to this script is 

~~~~~~~~~~~~~python
def init_my_case():
    n_compartment=500
    n_species=4
    liquid_phase = np.zeros((n_compartment,n_species))
    gas_phase = np.zeros((n_compartment,n_species)) # Not mandatory if only liquid 

    gas_phase[:,0]=0 # g/l
    gas_phase[:,1]=0.21 #g/l

    liquid_phase[:,0]=5 # g/l
    liquid_phase[:,1]=0 # g/l

    make_initial_concentration("./cma_data/my-case-init.h5",liquid_phase,gas_phase)
~~~~~~~~~~~~~


## CLI

Running the simulation directly through the command line with arguments passed to the executable is not recommended due to the complexity involved. Instead, the preferred approach is to use the runner Python script and provide a specific case descriptor file.
Below is a template for the case descriptor file, which should be in XML format:
~~~~~~~~~~~~~xml
<?xml version="1.0" encoding="UTF-8"?>
<!-- One file can contain multiple cases within this tag -->
<cases>
    <!-- 'name' refers to the case name to be given to the runner script -->
    <control name="my-case">
        <!-- Path to the root of the CMA case (don’t forget the leading '/') -->
        <!-- The 'recursive' parameter must be used for multi-flow map cases -->
        <cma_case_path recursive="true">/path/to/cma_case/</cma_case_path>
        
        <!-- Simulation time in seconds -->
        <final_time>10</final_time>
        
        <!-- Total number of particles to simulate (will be automatically balanced if using MPI) -->
        <number_particle>50000</number_particle>
        
        <!-- Simulation time step in seconds-->
        <!-- If 'n_compartment' > 1, leave the value as 0 to automatically determine the step size  -->
        <delta_time>0</delta_time>
        
        <!-- Name of the results folder -->
        <results_file_name>debug</results_file_name>
        
        <!-- Number of results to be saved in the HDF5 file -->
        <number_exported_result>100</number_exported_result>
        
        <!-- Model name; if not specified, the default model will be used -->
        <model_name>_</model_name>
        
        <!-- Scalar initializer file -->
        <initialiser_path>./cma_data/my-case-init.h5</initialiser_path>
    </control>
</cases>
~~~~~~~~~~~~~