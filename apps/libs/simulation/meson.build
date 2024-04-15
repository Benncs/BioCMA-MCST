
src_dir = './src'

simu_sources = run_command('sh','-c', f'find @src_dir@ -type f -name *.cpp ').stdout().splitlines()

public_simu  =  include_directories('public')



public_include_dir = [public_simu]

project_dependencies = [biomodel_dependency,mc_dependency,common_dependecy]
external_dependcies = [omp,cma_dep,blas,eigen_dep]
dep = [project_dependencies,external_dependcies]

simu_lib = library(f'@project_name@_simulation', simu_sources,
include_directories :[public_include_dir, 'includes'],dependencies:dep,pic:true)

simu_dependency = declare_dependency(
   include_directories : [public_include_dir],
   link_with : simu_lib,dependencies:dep)


subdir('tests')