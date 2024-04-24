

   src_dir = './src'
   message('Dynamic module feature activated')
   py = import('python').find_installation(pure: false)
   pybind11_dep = dependency('pybind11', required: true)
   pydep = py.dependency(embed:true)

   python_module_name = 'pyBioCMAMCST'
   modfile = files(f'@src_dir@/pymodule.cpp',f'@src_dir@/opaque_type.cpp')

   python_sources = files(f'@root@/modules/simple_model.py')

   # Install Python sources to 'birem' subdirectory
   py.install_sources(python_sources, subdir: python_module_name)

   install_subdir(
   'public',
   install_dir: includes_install_dir,
)


   dependencies = [
      biomodel_dependency,
      cma_dep,
      mc_dependency,
      common_dependecy,
      thread_dep,
      pydep,pybind11_dep,
   ]
   public = './public'

   indir = [public]

   py.extension_module(
      python_module_name,
      modfile,
      include_directories: indir,
      dependencies: [dependencies],
      install: true,subdir: python_module_name
   )

   py_dynamic_module = declare_dependency(
      include_directories: indir,
      dependencies: [pybind11_dep, pydep],
      sources:[f'@src_dir@/import_py.cpp'],compile_args:[use_python_define]
   )

