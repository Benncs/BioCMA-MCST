project(
  'biocma_mcst',
  'cpp',
  version: '0.1.1',
  meson_version: '>= 1.3.0',
  default_options: ['warning_level=3', 'b_colorout=always', 'cpp_std=gnu++20'],
)

use_cuda = get_option('use_cuda')

arch_specific_args = []

compiler = meson.get_compiler('cpp')

if compiler.has_argument('-fno-math-errno')
  arch_specific_args = [arch_specific_args, '-fno-math-errno']
endif

flto_flag = ''

if compiler.has_argument('-flto')
  flto_flag = '-flto'
endif

if compiler.has_argument('-mtune=native')
  arch_specific_args = [arch_specific_args, '-mtune=native']
endif

if compiler.has_argument('-Wno-interference-size')
  # Disable cache line detection size 
  arch_specific_args = [arch_specific_args, '-Wno-interference-size']
endif

compiler.has_argument('-fopenmp')

## Define main shared variables

outdir = meson.current_build_dir()
root = meson.project_source_root()
project_name = meson.project_name()
build_type = get_option('buildtype')
cpp = meson.get_compiler('cpp')
global_install_dir = get_option('prefix') + '/' + get_option('libdir')
includes_install_dir = ('include/' + project_name)

cpp_args = [arch_specific_args,'-DENABLE_KOKKOS_PROFILING']

add_project_arguments(cpp_args, language: 'cpp')

## Get buildtype specific options

CMAKE_TYPE='Release'
if build_type == 'debug' or build_type == 'debugoptimized'
  # NOTHING
  CMAKE_TYPE='Debug'
  add_global_arguments('-DDEBUG', language: 'cpp')
else
  add_global_arguments('-DNDEBUG', language: 'cpp')
  if compiler.has_argument('-ffast-math')
    arch_specific_args = [arch_specific_args, '-ffast-math']
  endif
endif

no_mpi_flag = '-DNO_MPI'
no_use_mpi_flag = '-DNO_USE_MPI'

use_dynamic_module = get_option('use_dynamic_module') # Use python module or not


use_python_define = '-DUSE_PYTHON_MODULE' # Define flag to pass to compiler if  use_dynamic_module is true


## Checking all needed depencies (execpt Python)


mpi = dependency('mpi', language: 'cpp', required: true)
thread_dep = dependency('threads', required: true)
omp = dependency('openmp')
hwloc = dependency('hwloc')
eigen_dep = dependency('eigen3', required: true)
blas = dependency('blas', required: false) # Eigen can work without BLAS


# if blas.found()
#   add_global_arguments('-DEIGEN_USE_BLAS', language: 'cpp')
# endif

#add_global_arguments('-DEIGEN_DONT_PARALLELIZE',language:'cpp')

__cmake_kokkos_dep__ = dependency(
  'Kokkos',
  method: 'cmake',
  required: true,
  components: ['kokkos'],
  modules: ['Kokkos::kokkos'],
)

kokkos_compile_arg = ['-DKOKKOS_DEPENDENCE']

if use_cuda
  kokkos_compile_arg = [
    kokkos_compile_arg,
    '-extended-lambda',
    '-Wext-lambda-captures-this',
    '-arch=native',
  ]
endif

kokkos = declare_dependency(
  dependencies: [__cmake_kokkos_dep__, omp, hwloc, thread_dep],
  compile_args: kokkos_compile_arg,
)

dynlib_dep = dependency('dynlib-1.0', required: true)

highfive_found = false
ceral_found = false
cmake = import('cmake', required: false)
if cmake.found()
  highfive_cmake_options = cmake.subproject_options()
  highfive_cmake_options.add_cmake_defines({'CMAKE_CXX_COMPILER': 'clang++-18'})
  highfive_cmake_options.add_cmake_defines({'CMAKE_CXX_STANDARD': '20'})

  highfive_cmake_options.add_cmake_defines({'CMAKE_BUILD_TYPE': CMAKE_TYPE})
  highfive_cmake_options.add_cmake_defines({'HIGHFIVE_USE_BOOST': false})
  highfive_cmake_options.add_cmake_defines({'HIGHFIVE_EXAMPLES': false})
  highfive_cmake_options.add_cmake_defines({'HIGHFIVE_UNIT_TESTS': false})
  highfive_cmake_options.add_cmake_defines({'HIGHFIVE_USE_BOOST': false})
  ## MAYBE USE ?? IDK 
  highfive_cmake_options.add_cmake_defines({'HDF5_PREFER_PARALLEL': false})
  highfive_cmake_options.add_cmake_defines({'HDF5_IS_PARALLEL': false})
  highfive_proj = cmake.subproject('highfive', options: highfive_cmake_options, required: false)
  if highfive_proj.found()
    highfive_found = true
    add_project_arguments('-DUSE_HIGHFIVE', language: 'cpp')
    highfive_dependency = highfive_proj.dependency('HighFive')
  endif

  cereal_cmake_options = cmake.subproject_options()
  cereal_cmake_options.add_cmake_defines({'SKIP_PERFORMANCE_COMPARISON': true})
  cereal_cmake_options.add_cmake_defines({'CMAKE_BUILD_TYPE': CMAKE_TYPE})
  cereal_cmake_options.add_cmake_defines({'BUILD_DOC': false})
  cereal_cmake_options.add_cmake_defines({'BUILD_SANDBOX': false})
  cereal_cmake_options.add_cmake_defines({'BUILD_SANDBOX': false})
  cereal_cmake_options.add_cmake_defines({'BUILD_TESTS': false})
  cereal_proj = cmake.subproject('cereal', options: cereal_cmake_options, required: false)
  if cereal_proj.found()
  ceral_found=true
  add_project_arguments('-DUSE_CEAREAL', language: 'cpp')
  cereal_dependency = cereal_proj.dependency('cereal')
  endif 

endif

cmt = dependency('cmtool', required: false)
f_cmt = false

if cmt.found()
  f_cmt = true
  cma_read_name = cmt.get_variable('cma_read', default_value: 'none')

  cma_dep = dependency(cma_read_name, required: false)

  if not cma_dep.found()
    f_cmt = false
  endif
endif

if f_cmt == false
  # Retrieve CMTOOL source code

  # TODO: remove core build 
  cma_reader = subproject(
    'cmtool',
    required: true,
    default_options: [
      'warning_level=3',
      'b_colorout=always',
      'b_coverage=false',
      'build_all=false',
      'python_wrap=false',
      'build_executable=false',
      'build_test=false',
    ],
  ) # We only need native cmtool library in order to read flowmap so, no python neither executable


  cma_dep = cma_reader.get_variable('cma_read_shared_lib_dependency') # Get shared lib


endif
subdir('apps')


run_command('sh', '-c', f'ln -s @outdir@/compile_commands.json .', check: false)
