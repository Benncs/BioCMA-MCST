src_dir = './src'

lib_sources = run_command('sh', '-c', f'find @src_dir@ -type f -name *.cpp ').stdout().splitlines()

name = f'@project_name@_api'

api_source = files(f'@src_dir@/api.cpp')

lib_api_shared = static_library(
  f'@name@_shared',
  api_source,
  include_directories: ['./includes', 'public'],
  dependencies: [core_shared_dependency],
  # cpp_args: ['-DNO_MPI'],
  install_rpath: project_name,
)

lib_api_distributed = static_library(
  f'@name@',
  api_source,
  include_directories: ['./includes', 'public'],
  dependencies: [core_distributed_dependency_static,mpi_wrap_dependency],
  # cpp_args: ['-DNO_MPI'],
  install_rpath: project_name,
)

lib_api_shared_dep = declare_dependency(
  link_with: lib_api_shared,
  include_directories: include_directories('./includes', 'public'),
  dependencies: [core_shared_dependency,kokkos],
)

lib_api_distributed_dep = declare_dependency(
  link_with: lib_api_distributed,
  include_directories: include_directories('./includes', 'public'),
  dependencies: [core_distributed_dependency_static],
)

build_python = get_option('build_python_api')
build_c_api = get_option('build_c_api')
if build_python
build_c_api = true # python API depends on C API 
endif 

# C API

if build_c_api

lib_api_c = both_libraries(
  f'@name@_c',
  files(f'@src_dir@/api_raw.cpp'),
  include_directories: ['./includes', 'public'],
  dependencies: [lib_api_distributed_dep],
  install_rpath: project_name,
)

lib_api_c_shared=library(
  f'@name@_c_shared',
  files(f'@src_dir@/api_raw.cpp'),
  include_directories: ['./includes', 'public'],
  dependencies: [lib_api_shared_dep],
  install_rpath: project_name,
)

lib_api_c_dep = declare_dependency(
  link_with: lib_api_c,
  include_directories: include_directories('./includes', 'public'),
)

lib_api_c_dep_static = declare_dependency(
  link_with: lib_api_c.get_static_lib(),
  include_directories: include_directories('./includes', 'public'),
)

lib_api_c_shared_dep = declare_dependency(
  link_with: lib_api_c_shared,
  include_directories: include_directories('./includes', 'public'),
)

# Python 
# Depends on C API 

if build_python


py = import('python').find_installation(pure: false)
pybind11_dep = dependency('pybind11', required: true)
pydep = py.dependency(embed: true)
python_module_name = 'handle_module'
modfile = files(f'@src_dir@/api_python.cpp')

py.extension_module(
  python_module_name,
  modfile,
  include_directories: ['./includes', 'public'],
  install: true,
  dependencies: [lib_api_distributed_dep, lib_api_c_dep_static, pydep, pybind11_dep],
  subdir: python_module_name,
)

python_sources = files(f'./python/__init__.py')
py.install_sources(python_sources, subdir: python_module_name)
endif 

endif 

if build_test
subdir('tests')
endif
