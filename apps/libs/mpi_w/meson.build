# This sublibrary aims to provide wrap around MPI C API for modern C++


src_dir = './src'

mpi_wrap_sources = run_command('sh', '-c', f'find @src_dir@ -type f -name *.cpp ').stdout().splitlines()

public_mpi_wrap = include_directories('public')

public_mpi_w_include_dir = [public_mpi_wrap]
dep = [mpi,common_dependecy,cma_dep]

mpi_wrap_lib = library(
   f'@project_name@_mpi_wrap',
   mpi_wrap_sources,
   include_directories: [public_mpi_w_include_dir],
   dependencies: dep,
   pic: true,
   install: true,
   install_rpath: global_install_dir,
)

mpi_wrap_dependency = declare_dependency(
   include_directories: public_mpi_w_include_dir,
   dependencies:dep,
   link_with: mpi_wrap_lib,
)

install_subdir(
   'public',
   install_dir: includes_install_dir,
)

subdir('./tests')