

src_dir = './src'

app_sources = run_command('sh', '-c', f'find @src_dir@ -type f -name *.cpp ').stdout().splitlines()


lib_sources = []
foreach item : app_sources
  if item not in [f'@src_dir@/main.cpp', f'@src_dir@/import_py.cpp']
    lib_sources += item
  endif
endforeach




dep_all = [
  mpi_wrap_dependency,
  biomodel_dependency,
  cma_dep,
  mc_dependency,
  common_dependecy,
  thread_dep,
]




if highfive_found
  dep_all += highfive_proj.dependency('HighFive')
endif



dep_static_module = [
  dep_all,
  simulation_lib_dependency,
  omp,
]



name = f'@project_name@_cli_app'


# Define static library and its dependency
lib_app = static_library(
  name,
  lib_sources,
  include_directories: ['./includes'],
  dependencies: [dep_static_module],
  cpp_args: [],
  install_rpath: project_name,
)

app_test_dependency = declare_dependency(
  include_directories: ['./includes'],
  link_with: lib_app,
  dependencies: dep_static_module,
  compile_args: [],
)

# Define executable
executable(
  name,
  f'@src_dir@/main.cpp',
  include_directories: ['./includes'],
  dependencies: [app_test_dependency],
  cpp_args: [target_custom_options],
  install: true,
  install_rpath: global_install_dir,
  link_args: [flto_flag]
)

# Conditional dynamic module build
if use_dynamic_module
  dep_dynamic_module = [dep_all, simu_dyn_dependency, py_dynamic_module]
  name = name + '_dyn'

  lib_app_dyn = static_library(
    name,
    lib_sources,
    include_directories: ['./includes'],
    dependencies: dep_dynamic_module,
    cpp_args: [use_python_define],
  )

  app_test_dependency_dyn = declare_dependency(
    include_directories: ['./includes'],
    link_with: lib_app_dyn,
    dependencies: dep_dynamic_module,
    compile_args: [use_python_define],
  )

  executable(
    name,
    f'@src_dir@/main.cpp',
    include_directories: ['./includes'],
    dependencies: app_test_dependency_dyn,
    cpp_args: [use_python_define],
    install: true,
  )
endif


subdir('tests')

