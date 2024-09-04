# Directory containing simulation source files
simulation_src_dir = './src'

# Find all source files in the simulation source directory
simulation_sources = run_command(
   'sh',
   '-c',
   f'find @simulation_src_dir@ -type f -name *.cpp',
).stdout().splitlines()

lib_name = f'@project_name@_simulation'

public_simulation_includes = include_directories('public')

install_subdir(
   'public',
   install_dir: includes_install_dir,
)



private_simulation_includes = include_directories('includes')

project_dependencies = [biomodel_dependency, mc_dependency, common_dependecy]

# Private dependencies needed to compile .so 
private_dependencies = [eigen_dep]

external_dependcies = [cma_dep, blas,kokkos]

simulation_public_dependencies = [
   project_dependencies,
   external_dependcies,
]

simulation_dyn_public_dep = [project_dependencies, external_dependcies]

# Create the simulation library
simulation_lib = static_library(
   lib_name,
   simulation_sources,
   include_directories: [
      public_simulation_includes,
      private_simulation_includes,
   ],
   dependencies: [simulation_public_dependencies, private_dependencies],
   # pic: true,
   install: true,
   install_rpath: global_install_dir,
   link_args: [flto_flag],
   cpp_args:[]
)


simulation_lib_dependency = declare_dependency(
   include_directories: [public_simulation_includes],
   link_with: simulation_lib,
   dependencies: simulation_public_dependencies,
)



# Check if the dynamic module is enabled
if use_dynamic_module
   # Create the dynamic simulation library
   simulation_lib_dyn = static_library(
      f'@lib_name@_dyn',
      simulation_sources,
      include_directories: [
         public_simulation_includes,
         private_simulation_includes,
      ],
      dependencies: [simulation_dyn_public_dep, private_dependencies],
      cpp_args: [use_python_define],
      # pic: true,
      install: true,
       install_rpath: global_install_dir
   )

   # Declare dependency for the dynamic simulation library
   simu_dyn_dependency = declare_dependency(
      include_directories: [public_simulation_includes],
      link_with: simulation_lib_dyn,
      dependencies: simulation_dyn_public_dep,
      compile_args: [use_python_define],
      
   )
endif



subdir('tests')
