cma_utils_src = run_command(
   'sh',
   '-c',
   f'find ./src -type f -name *.cpp',check:true
).stdout().splitlines()



src_without_cmtool_dep  = run_command(
   'sh',
   '-c',
   f'find ./transitionner -type f -name *.cpp',check:true
).stdout().splitlines()



lib_name = f'@project_name@_cma_utils'

public = include_directories('public')
private = include_directories('includes')


deps = [eigen_dep,blas, blas,kokkos,common_dependecy]


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

include_transitionner = include_directories('transitionner')
transitionner_lib = static_library(
   f'@project_name@_transitionner',
   src_without_cmtool_dep,
   include_directories: [
      public,private,include_transitionner
   ],
   dependencies: [deps,cma_dep],
   install_rpath: global_install_dir,
   link_args: [flto_flag],
   cpp_args:[]
)



transitionner_dependency = declare_dependency(
   include_directories: [public,include_transitionner],
   link_with: transitionner_lib,
   dependencies: [deps,cma_dep],
)


cma_utils_dependency = declare_dependency(
   include_directories: [public],
   link_with: cma_utils,
   dependencies: deps,
)
