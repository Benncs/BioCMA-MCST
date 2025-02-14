cma_utils_src = run_command(
   'sh',
   '-c',
   f'find ./src -type f -name *.cpp',
).stdout().splitlines()

lib_name = f'@project_name@_cma_utils'

public = include_directories('public')
private = include_directories('includes')


deps = [eigen_dep,blas,cma_dep, blas,kokkos,common_dependecy]



cma_utils = static_library(
   lib_name,
   cma_utils_src,
   include_directories: [
      public,private
   ],
   dependencies: deps,
   # pic: true,
   # install: true,
   install_rpath: global_install_dir,
   link_args: [flto_flag],
   cpp_args:[]
)


cma_utils_dependency = declare_dependency(
   include_directories: [public],
   link_with: cma_utils,
   dependencies: deps,
)
