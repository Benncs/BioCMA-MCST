src_dir = './src'
message('Dynamic module feature activated')
py = import('python').find_installation(pure: false)
pybind11_dep = dependency('pybind11', required: true)
pydep = py.dependency(embed: true)

python_module_name = 'pyBioCMAMCST'
modfile = files(f'@src_dir@/pymodule.cpp', f'@src_dir@/opaque_type.cpp')

# python_sources = files(f'@root@/modules/simple_model.py')

# Install Python sources to 'birem' subdirectory
# py.install_sources(python_sources, subdir: python_module_name)

install_subdir(
  'public',
  install_dir: includes_install_dir,
)
public = include_directories('./public')

# dependencies = [
#   biomodel_dependency,
#   cma_dep,
#   mc_dependency,
#   common_dependecy,
#   pydep,
#   pybind11_dep,
# ]

# py.extension_module(
#   python_module_name,
#   modfile,
#   include_directories: public,
#   dependencies: [dependencies],
#   install: true,
#   subdir: python_module_name,
# )

wrap_py = library(
  'wrap_py',
  include_directories: public,
  sources: [f'@src_dir@/import_py.cpp',f'@src_dir@/pymodule.cpp'],
  dependencies: [pybind11_dep, pydep, mc_dependency],
)

py_dynamic_module = declare_dependency(include_directories: public, link_with: wrap_py)
