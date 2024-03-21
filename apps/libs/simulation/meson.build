
src_dir = './src'

simu_sources = run_command('sh','-c', f'find @src_dir@ -type f -name *.cpp ').stdout().splitlines()

public_simu  =  include_directories('public')

public_include_dir = [public_simu,public_common,eigen_includes,]

simu_lib = static_library(f'@project_name@_simulation', simu_sources,
include_directories :[public_include_dir, 'includes'],dependencies:[omp,mc_dependency])

simu_dependency = declare_dependency(
   include_directories : [public_include_dir,'includes'],
   link_with : simu_lib,
   sources : simu_sources)


