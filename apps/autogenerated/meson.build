# Use to generate source/header file before compilation 

version = meson.project_version().split('.')
assert(version.length()==3)

conf_data = configuration_data(
    {
        '_BIOMC_VERSION_MAJOR': version[0],
        '_BIOMC_VERSION_MINOR': version[1],
        '_BIOMC_VERSION_DEV':  version[2],
        'probe_buffer_size': '10e3',
        '__f__redirect__': 'false',
        '__f_verbose__': 'true',
        'dump_particle_state': 'true',
        'use_probe': 'true',
        'const_number_simulation': 'false',
        'default_particle_container_allocation_factor': '2.5',
        'MC_MAX_PARTICLE_BUFFER': '100e6',
        'debug_MC_RAND_DEFAULT_SEED': '2024',
        '_dead_particle_ratio_threshold':'0.01',
        'result_dir_path':'"./results/"',
        'cma_data_folder_path':'"cma_data"'
    },
)

configure_file(
    input: 'biocma_cst_config.hpp.in',
    output: 'biocma_cst_config.hpp',
    configuration: conf_data,
)

prog_python = find_program('python3')

model_loader = meson.current_source_dir()+'/model_loader.cxx'
model_loader_header = meson.current_source_dir()+'/model_loader.hxx'
variant_model = meson.current_source_dir()+'/variant_model.hxx'

arg_includ_py=''
if use_dynamic_module
arg_includ_py='1'
endif

arg_udf='1'
if use_cuda
arg_udf=''
endif


# custom target to generate model_loader.cpp and model_list.hpp
generated_model_selection = custom_target(
    'generated_model_selection.[ch]',
    output: [
        'wrap_init_model_selector.cpp',
        'wrap_init_model_selector.hpp',
        'variant_model.hpp',
    ],
    input: f'@root@/devutils/auto_generate/generate_model_loader.py',
    command: [
        prog_python,
        '@INPUT@',
        f'@root@/apps/libs/biomodels',
        model_loader,
        '@OUTPUT0@',
        model_loader_header,
        '@OUTPUT1@',
        variant_model,
        '@OUTPUT2@',
        arg_includ_py,arg_udf
    ],
)