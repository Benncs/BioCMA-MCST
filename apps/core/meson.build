src_dir = './src'

core_source = run_command('sh', '-c', f'find @src_dir@ -type f -name *.cpp ').stdout().splitlines()

dep_core = [
    kokkos,
    biomodel_dependency,
    mc_dependency,
    common_dependecy,
    model_selection_dep,
    simulation_lib_dependency,
    eigen_dep,blas,
    transitionner_dependency
]

if highfive_found
    dep_core += highfive_dependency
endif

if ceral_found
    dep_core += cereal_dependency
endif

dep_core_mpi = [dep_core,mpi_wrap_dependency]


public_core_includes = include_directories('public')

private_core_includes = include_directories('includes')

core_shared = static_library(
    'core_shared',
    core_source,
    include_directories: [private_core_includes, public_core_includes],
    dependencies: [dep_core],
    cpp_args: [no_mpi_flag],
    # pic: true,
)

core_distributed = static_library(
    'core_distributed',
    core_source,
    include_directories: [private_core_includes, public_core_includes],
    dependencies: [dep_core_mpi], 
    # pic: true
)

core_shared_dependency = declare_dependency(
    include_directories: [public_core_includes],
    link_with: core_shared,
    dependencies: dep_core,
    compile_args: [no_mpi_flag],
)

core_distributed_dependency = declare_dependency(
    include_directories: [public_core_includes],
    link_with: core_distributed,
    dependencies: dep_core_mpi,
    compile_args: [],
)

core_distributed_dependency_static = declare_dependency(
    include_directories: [public_core_includes],
    link_with: core_distributed,
    dependencies: dep_core_mpi,
    compile_args: [],
)

subdir('tests')