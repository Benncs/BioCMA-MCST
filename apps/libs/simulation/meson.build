
src_dir = './src'

simu_sources = run_command('sh','-c', f'find @src_dir@ -type f -name *.cpp ').stdout().splitlines()

public_simu  =  include_directories('public')


private_simu_includes = include_directories('./includes')
public_include_dir = [public_simu]

project_dependencies = [biomodel_dependency,mc_dependency,common_dependecy]
private_dep = [eigen_dep]
external_dependcies = [omp,cma_dep,blas,eigen_dep]
simulation_public_dep = [project_dependencies,external_dependcies]

simu_lib = library(f'@project_name@_simulation', simu_sources,
include_directories :[public_include_dir, private_simu_includes],dependencies:[simulation_public_dep,private_dep],pic:true)

simu_dependency = declare_dependency(
   include_directories : [public_include_dir],
   link_with : simu_lib,dependencies:simulation_public_dep)


subdir('tests')