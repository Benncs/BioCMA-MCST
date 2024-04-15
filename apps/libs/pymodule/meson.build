src_dir = './src'
use_dynamic_module = get_option('use_dynamic_module')

if use_dynamic_module
   python_extension_suffix = run_command('python3-config', '--extension-suffix').stdout().strip()
   python_extension_suffix = python_extension_suffix.substring(1)
   python_module_name = 'pyBioCMAMCST'
   modfile = files(f'@src_dir@/pymodule.cpp') 


   dependencies = [biomodel_dependency,cma_dep,mc_dependency,common_dependecy,thread_dep,python_embded]
   public='./public'

   indir = [public,'/usr/include/python3.10/']

   py = shared_library(python_module_name,
                     modfile,
                     name_suffix: python_extension_suffix,
                     name_prefix: '',
                     include_directories: indir,
                     dependencies: [dependencies],link_args : ['-L/usr/lib', '-lpython3.10'])

   py_dynamic_module = declare_dependency(
      include_directories : indir,dependencies:python_embded,link_args : ['-L/usr/lib', '-lpython3.10'])
endif 
