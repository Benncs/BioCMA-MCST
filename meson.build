project(
  'biocma_mcst',
  'cpp',
  version: '0.0.4',
  meson_version: '>= 1.3.0',
  default_options: [
    'warning_level=3',
    'b_colorout=always',
  ],
)

use_cuda = get_option('use_cuda')

arch_specific_args = []
if use_cuda
  message('USING CUDA')
  arch_specific_args = ['-mp=gpu' ,'-gpu=cc35','-cudalib=cublas']
  # add_project_arguments('cpp_std=gnu++20',language: 'cpp')
  
else
  arch_specific_args = ['-march=native']
  # arch_specific_args=['-mavx2']
  cpp_version = get_option('cpp_std')
  if cpp_version != 'gnu++20'
    error('c++20 features should be available')
  endif 
endif 
compiler = meson.get_compiler('cpp')

if compiler.has_argument('-fno-math-errno')
arch_specific_args = [arch_specific_args,'-fno-math-errno']
endif 

flto_flag = ''

if compiler.has_argument('-flto')
flto_flag = '-flto'
endif 

if compiler.has_argument('-mtune=native')
arch_specific_args = [arch_specific_args,'-mtune=native']
endif 





compiler.has_argument('-fopenmp')

## Define main shared variables\

compiler.has_argument('--std=c++20',required: true)
add_global_arguments('--std=c++20',language: 'cpp')

root = meson.project_source_root()
project_name = meson.project_name()
build_type = get_option('buildtype')
cpp = meson.get_compiler('cpp')

# global_install_dir = 'lib' # 'lib/'+project_name  # Relative  to /usr/local/
global_install_dir = get_option('prefix')+'/'+get_option('libdir')
includes_install_dir = ('include/' + project_name)

cpp_args = arch_specific_args 


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
  # if compiler.has_argument('-pg')
  # arch_specific_args= [arch_specific_args,'-pg']
  # add_project_link_arguments('-pg',language: 'cpp')
  # endif 

else
  add_global_arguments('-DNDEBUG', language: 'cpp')
    if compiler.has_argument('-ffast-math')
  arch_specific_args= [arch_specific_args,'-ffast-math']
  endif 
endif

use_dynamic_module = get_option('use_dynamic_module') # Use python module or not
use_python_define = '-DUSE_PYTHON_MODULE' # Define flag to pass to compiler if  use_dynamic_module is true

## Checking all needed depencies (execpt Python)

mpi = dependency('mpi', language: 'cpp', required: true)
thread_dep = dependency('threads', required: true)

if use_cuda
omp = declare_dependency(link_args: ['-fopenmp'])
else
omp = dependency('openmp', language: 'cpp', required: true)
endif 
tbb = dependency('tbb', required: true)


eigen_dep = dependency('eigen3', required: true)
blas = dependency('blas', required: false) # Eigen can work without BLAS


if blas.found()
  add_global_arguments('-DEIGEN_USE_BLAS', language: 'cpp')
endif

highfive_found = false
cmake = import('cmake', required: false)
if cmake.found()
  highfive_cmake_options = cmake.subproject_options()
  highfive_cmake_options.add_cmake_defines({'HIGHFIVE_EXAMPLES': false})
  highfive_cmake_options.add_cmake_defines({'HIGHFIVE_UNIT_TESTS': false})
  highfive_cmake_options.add_cmake_defines({'HDF5_IS_PARALLEL': true})
  highfive_cmake_options.add_cmake_defines({'CMAKE_CXX_STANDARD': '20'})
  highfive_cmake_options.add_cmake_defines({'HIGHFIVE_USE_BOOST': false})
  highfive_cmake_options.add_cmake_defines({'HDF5_PREFER_PARALLEL': true})
  highfive_proj = cmake.subproject('highfive', options: highfive_cmake_options, required: false)
  if highfive_proj.found()
    highfive_found = true
    add_project_arguments('-DUSE_HIGHFIVE', language: 'cpp')
    highfive_dependency = highfive_proj.dependency('HighFive')
  endif
endif

default_option_cpp20 = 'cpp_std=none'
if not use_cuda
  default_option_cpp20 = 'cpp_std=gnu++20'
endif 

cmt = dependency('cmtool', required: false)
f_cmt = false

if cmt.found()
  f_cmt = true
  cma_read_name = cmt.get_variable('cma_read',default_value :'none')

  cma_dep = dependency(cma_read_name, required: false)

  if not cma_dep.found()
  f_cmt=false
  endif
endif 

if f_cmt==false
# Retrieve CMTOOL source code
cma_reader = subproject(
  'cmtool',
  default_options: [
    'warning_level=3',
    'b_colorout=always',
    'b_coverage=false',
    'build_all=false',
    'python_wrap=false',
    'build_executable=false',
    'build_test=false',
    f'use_cuda=@use_cuda@',
    default_option_cpp20
  ],
) # We only need native cmtool library in order to read flowmap so, no python neither executable

cma_dep = cma_reader.get_variable('cma_read_shared_lib_dependency') # Get shared lib

endif 
subdir('apps')

# Clangd project configuration
run_command('sh', '-c', 'ln -s builddir/compile_commands.json .', check: false)
