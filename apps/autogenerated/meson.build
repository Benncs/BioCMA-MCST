# Use to generate source/header file before compilation

version = meson.project_version().split('.')
assert(version.length() == 3)

conf_data = configuration_data(
  {
    '_BIOMC_VERSION_MAJOR': version[0],
    '_BIOMC_VERSION_MINOR': version[1],
    '_BIOMC_VERSION_DEV': version[2],
    'probe_buffer_size': '10e3',
    '__f__redirect__': 'false',
    '__f_verbose__': 'true',
    'dump_particle_state': 'true',
    'use_probe': 'false',
    'const_number_simulation': 'false',
    'default_particle_container_allocation_factor': '2.5',
    'MC_MAX_PARTICLE_BUFFER': '100e6',
    'debug_MC_RAND_DEFAULT_SEED': '2024',
    '_dead_particle_ratio_threshold': '0.01',
    'cma_data_folder_path': '"cma_data"',
    'enable_event_counter': 'true',
  },
)
result_dir = get_option('result_dir_path')
conf_data.set_quoted('result_dir_path', result_dir)
conf_data.set_quoted('_BIOMC_BUILD_MODE', build_type)

configure_file(
  input: 'biocma_cst_config.hpp.in',
  output: 'biocma_cst_config.hpp',
  configuration: conf_data,
)

prog_python = find_program('python3')

template_model_loader = meson.current_source_dir() + '/model_loader.cxx'
template_model_loader_header = meson.current_source_dir() + '/model_loader.hxx'
template_variant_model = meson.current_source_dir() + '/variant_model.hxx'

arg_includ_py = ''
if use_dynamic_module
  arg_includ_py = '1'
endif

arg_udf = ''
if not use_cuda and build_udf
  arg_udf = '1'
endif

# model_list_name = ['two_meta','simple_model','two_meta_nb','fixed_length']
model_list_name = ['fixed_length','two_meta_nb']
if not use_cuda and build_udf
  model_list_name = [model_list_name, 'udf_model']
endif

model_folder_root = f'@root@/apps/libs/models'

generated_model_selection = custom_target(
  'generated_model_selection.[ch]',
  output: [
    'wrap_init_model_selector.cpp',
    'wrap_init_model_selector.hpp',
    'variant_model.hpp',
  ],
  input: f'@root@/devutils/generate_model_loader.py',
  command: [
    prog_python,
    '@INPUT@',
    model_folder_root,
    template_model_loader,
    template_model_loader_header,
    template_variant_model,
    '@OUTPUT0@',
    '@OUTPUT1@',
    '@OUTPUT2@',
    arg_includ_py,
    arg_udf,
    model_list_name,
  ],
)
