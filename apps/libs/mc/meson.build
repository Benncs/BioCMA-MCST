
src_dir = './src'

mc_sources = run_command('sh', '-c', f'find @src_dir@ -type f -name *.cpp ').stdout().splitlines()

public_mc = include_directories('public')

public_include_dir = [public_mc, public_common]

mc_lib = library(
   f'@project_name@_mc',
   mc_sources,
   include_directories: [public_include_dir, 'includes'],
   dependencies: [omp,cma_dep],
   pic: true,
   install: true,
   install_dir: global_install_dir,
)

mc_dependency = declare_dependency(
   include_directories: public_include_dir,
   dependencies:[omp,cma_dep],
   link_with: mc_lib,
)

install_subdir(
   'public',
   install_dir: includes_install_dir,
)

