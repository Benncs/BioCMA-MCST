v = meson.project_version().split('.')
assert(v.length()==3)

conf_data = configuration_data(
    {
        '_BIOMC_VERSION_MAJOR': v[0],
        '_BIOMC_VERSION_MINOR': v[1],
        '_BIOMC_VERSION_DEV':  v[2],
        'probe_buffer_size': '10e3',
        '__f__redirect__': 'false',
        '__f_verbose__': 'true',
        'dump_particle_state': 'true',
        'use_probe': 'true',
        'const_number_simulation': 'false',
        'default_particle_container_allocation_factor': '2.5',
        'MC_MAX_PARTICLE_BUFFER': '100e6',
        'debug_MC_RAND_DEFAULT_SEED': '2024',
    },
)

configure_file(
    input: 'biocma_cst_config.hpp.in',
    output: 'biocma_cst_config.hpp',
    configuration: conf_data,
)

prog_python = find_program('python3')

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
        f'@root@/devutils/auto_generate/model_loader.cxx',
        '@OUTPUT0@',
        f'@root@/devutils/auto_generate/model_loader.hxx',
        '@OUTPUT1@',
        f'@root@/devutils/auto_generate/variant_model.hxx',
        '@OUTPUT2@',
    ],
)