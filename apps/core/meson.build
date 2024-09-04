src_dir = './src'

core_source = run_command('sh', '-c', f'find @src_dir@ -type f -name *.cpp ').stdout().splitlines()

dep_core = [
    kokkos,
    biomodel_dependency,
    cma_dep,
    mc_dependency,
    common_dependecy,
    thread_dep,
    model_selection_dep,
    cereal_dependency,
    simulation_lib_dependency
]

if highfive_found
    dep_core += highfive_proj.dependency('HighFive')
endif


dep_core_mpi = [dep_core,mpi_wrap_dependency]


public_core_includes = include_directories('public')

core_shared = static_library(
    'core_shared',
    core_source,
    include_directories: ['./includes', public_core_includes],
    dependencies: [dep_core],
    cpp_args: ['-DNO_MPI'],
)

core_shared_dependency = declare_dependency(
    include_directories: [public_core_includes],
    link_with: core_shared,
    dependencies: dep_core,
    compile_args: ['-DNO_MPI'],
)

core_distributed = static_library(
    'core_distributed',
    core_source,
    include_directories: ['./includes', public_core_includes],
    dependencies: [dep_core_mpi],)

core_distributed_dependency = declare_dependency(
    include_directories: [public_core_includes],
    link_with: core_distributed,
    dependencies: dep_core_mpi,
    compile_args: [],
)