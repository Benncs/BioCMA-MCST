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
  dependencies: [core_distributed_dependency],
  # cpp_args: ['-DNO_MPI'],
  install_rpath: project_name,
)

lib_api_shared_dep = declare_dependency(
  link_with: lib_api_shared,
  include_directories: include_directories('./includes', 'public'),
  dependencies: [core_shared_dependency],
)

lib_api_distributed_dep = declare_dependency(
  link_with: lib_api_distributed,
  include_directories: include_directories('./includes', 'public'),
  dependencies: [core_distributed_dependency],
)

# C API

lib_api_c = library(
  f'@name@_c',
  files(f'@src_dir@/api_raw.cpp'),
  include_directories: ['./includes', 'public'],
  dependencies: [lib_api_distributed_dep],
  install_rpath: project_name,
)

library(
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

# Python 
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
  dependencies: [lib_api_distributed_dep, lib_api_c_dep, pydep, pybind11_dep],
  subdir: python_module_name,
  link_with: [lib_api_c, lib_api_distributed],
)

python_sources = files(f'@root@/modules/simple_model.py')
py.install_sources(python_sources, subdir: python_module_name)

install_subdir(
  'public',
  install_dir: includes_install_dir,
)
