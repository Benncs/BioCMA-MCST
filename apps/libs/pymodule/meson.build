src_dir = './src'
message('Dynamic module feature activated')
py = import('python').find_installation(pure: false)
pybind11_dep = dependency('pybind11', required: true)
pydep = py.dependency(embed: true)

install_subdir(
  'public',
  install_dir: includes_install_dir,
)
public = include_directories('./public')


wrap_py = library(
  'wrap_python',
  include_directories: public,
  sources: [f'@src_dir@/import_py.cpp',f'@src_dir@/pymodule.cpp'],
  dependencies: [pybind11_dep, pydep, mc_dependency],
)

py_dynamic_module = declare_dependency(include_directories: public, link_with: wrap_py)
