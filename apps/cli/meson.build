src_dir = './src'

app_sources = run_command('sh', '-c', f'find @src_dir@ -type f -name *.cpp ', check: true).stdout().splitlines()
lib_sources = []
foreach item : app_sources
  if item not in [f'@src_dir@/main.cpp', f'@src_dir@/import_py.cpp']
    lib_sources += item
  endif
endforeach

name = f'@project_name@_cli_app'

includes = ['./includes']

exe_link_args=[flto_flag]
if not use_cuda and build_udf
exe_link_args=[exe_link_args]
endif

lib_app = static_library(
  name,
  lib_sources,
  include_directories: includes,
  dependencies: [lib_api_distributed_dep],
  cpp_args: [],
  install_rpath: project_name,
)

lib_app_shared = static_library(
  f'@name@_shared',
  lib_sources,
  include_directories: includes,
  dependencies: [lib_api_shared_dep],
  # cpp_args: ['-DNO_MPI'],
  install_rpath: project_name,
)

app_dependency = declare_dependency(
  include_directories: includes,
  link_with: lib_app,
  dependencies: lib_api_distributed_dep,
  compile_args: [],
)

app_shared_dependency = declare_dependency(
  include_directories: includes,
  link_with: lib_app_shared,
  dependencies: lib_api_shared_dep,
  compile_args: [],
)

# Define executable
executable(
  name,
  f'@src_dir@/main.cpp',
  include_directories: includes,
  dependencies: [app_dependency, udf_loader],
  install: true,
  install_rpath: global_install_dir,
  install_dir: global_install_dir,
  link_args: exe_link_args,
)

executable(
  f'@name@_shared',
  f'@src_dir@/main.cpp',
  include_directories: includes,
  dependencies: [app_shared_dependency, udf_loader],

  install: true,
  install_rpath: global_install_dir,
  install_dir: global_install_dir,
  link_args: exe_link_args,
)

# Conditional dynamic module build
if use_dynamic_module
  # dep_dynamic_module = [app_dependency, py_dynamic_module]
  # name = name + '_dyn'

  # lib_app_dyn = static_library(
  #   name,
  #   lib_sources,
  #   include_directories: includes,
  #   dependencies: dep_dynamic_module,
  #   cpp_args: [use_python_define],
  # )

  # app_test_dependency_dyn = declare_dependency(
  #   include_directories: includes,
  #   link_with: lib_app_dyn,
  #   dependencies: dep_dynamic_module,
  #   compile_args: [use_python_define],
  # )

  # executable(
  #   name,
  #   f'@src_dir@/main.cpp',
  #   include_directories: includes,
  #   dependencies: [app_test_dependency_dyn,udf_loader],
  #   cpp_args: [use_python_define],
  #   install: true,
  # )
endif

if build_test
  subdir('tests')
endif
