
src_dir = './src'

mc_sources = run_command('sh', '-c', f'find @src_dir@ -type f -name *.cpp ').stdout().splitlines()

public_mc = include_directories('public')

models_list_def = include_directories('../biomodels/public')

public_include_dir = [public_mc, public_common,models_list_def,'../']

mc_lib = static_library(
   f'@project_name@_mc',
   mc_sources,
   include_directories: [public_include_dir, 'includes','../'],
   dependencies: [cma_dep,kokkos],
   pic: true,
   install: true,
   install_rpath: global_install_dir,
   link_args: [flto_flag]
)

mc_dependency = declare_dependency(
   include_directories: public_include_dir,
   dependencies:[cma_dep,kokkos],
   link_with: mc_lib,
)

subdir('./tests')

install_subdir(
   'public',
   install_dir: includes_install_dir,
)

