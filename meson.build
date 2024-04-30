project(
  'biocma_mcst',
  'cpp',
  version: '0.0.2',
  meson_version: '>= 1.3.0',
  default_options: [
    'warning_level=3',
    'cpp_std=gnu++20',
    'b_colorout=always',
    # 'b_coverage=true'
  ],
)

## Define main shared variables 
root = meson.project_source_root()
project_name = meson.project_name()
build_type = get_option('buildtype')
cpp = meson.get_compiler('cpp')

global_install_dir = 'lib' # 'lib/'+project_name  # Relative  to /usr/local/
includes_install_dir = 'include/'+project_name # Relative  to /usr/local/

cpp_args = ['-march=native', '-fno-math-errno']
add_project_arguments(cpp_args, language: 'cpp')

## Get buildtype specific options 

#TODO FIXME 
project_custom_options = {
  'debug': [],
  'release': [],
}
target_custom_options = project_custom_options.get(build_type, [])
# foreach opt: target_custom_options
#     add_project_arguments(opt, language: 'cpp')
#     message(opt)
# endforeach

if build_type == 'debug' or build_type == 'debugoptimized'
  add_global_arguments('-DDEBUG', language: 'cpp')
endif



use_dynamic_module = get_option('use_dynamic_module') # Use python module or not 
use_python_define = '-DUSE_PYTHON_MODULE' # Define flag to pass to compiler if  use_dynamic_module is true 


## Checking all needed depencies (execpt Python) 

mpi = dependency('mpi', language: 'cpp', required: true)
thread_dep = dependency('threads', required: true)
omp = dependency('openmp', language: 'cpp', required: true)
eigen_dep = dependency('eigen3', required: true)

blas = dependency('blas', required: false) # Eigen can work without BLAS 
if blas.found()
  add_global_arguments('-DEIGEN_USE_BLAS', language: 'cpp')
endif


# Retrieve BIREM source code  
cma_reader = subproject(
  'birem',
  default_options: [
    'warning_level=3',
    'cpp_std=gnu++20',
    'b_colorout=always',
    'build_all=false',
    'python_wrap=false',
    'build_executable=false',
    'build_test=false',
  ]
) # We only need native birem library in order to read flowmap so, no python neither executable 

cma_dep = cma_reader.get_variable('cma_read_shared_lib_dependency') # Get shared lib 

subdir('apps')

# Clangd project configuration 
run_command('sh', '-c', 'ln -s builddir/compile_commands.json .', check: false)




