project(
    'biocma_mcst',
    'cpp',
    version: '1.2.1',
    meson_version: '>= 1.8.0',
    default_options: ['warning_level=3', 'b_colorout=always', 'cpp_std=c++20', 'b_pch=false'],
)

use_cuda = get_option('use_cuda')

if use_cuda and get_option('b_pch')
    error(
        'Precompiled header not supported for CUDA build\r\nPlease reconfigure the project turning b_pch to false\r\nmeson setup builddir --reconfigure -Db_pch=false',
    )
endif

build_type = get_option('buildtype')
python_install = get_option('python_install')
ci_execution = get_option('ci_execution')
use_dynamic_module = get_option('use_dynamic_module') # Use python module or not
use_cereal = get_option('use_cereal')
use_system_kokkos = get_option('use_system_kokkos')
build_test = get_option('build_test')
build_udf = get_option('build_udf')

if not use_cuda and build_udf
    add_project_arguments('-DUSE_UDF', language: 'cpp')
endif

if use_cuda and build_udf
    error(
        'UDF model is not supported for CUDA build\r\nPlease reconfigure the project turning build_udf to false\r\nmeson setup builddir --reconfigure -Dbuild_udf=false',
    )
endif

no_mpi_flag = '-DNO_MPI'

arch_specific_args = []

compiler = meson.get_compiler('cpp')

flto_flag = []

if compiler.has_argument('-flto')
    flto_flag = [flto_flag, '-flto']
endif

if compiler.has_argument('-dlto')
    flto_flag = [flto_flag, '-dlto']
endif

if compiler.has_argument('-march=native')
    arch_specific_args = [arch_specific_args, '-march=native']
endif

if compiler.has_argument('-mtune=native')
    arch_specific_args = [arch_specific_args, '-mtune=native']
endif

if compiler.has_argument('-Wno-interference-size')
    # Disable cache line detection size
    arch_specific_args = [arch_specific_args, '-Wno-interference-size']
endif

compiler.has_argument('-fopenmp')

#if (compiler.get_id() == 'clang')
#message('using clang')
## add_global_arguments(['-stdlib=libc++'], language: 'cpp')
## add_global_link_arguments(['-stdlib=libc++', '-lc++', '-lc++abi'], language: 'cpp')
#endif

## Define main shared variables
outdir = meson.current_build_dir()
root = meson.project_source_root()
project_name = meson.project_name()
cpp = meson.get_compiler('cpp')
if python_install
    global_install_dir = get_option('prefix') + '/' + get_option('libdir')
else
    global_install_dir = '/opt/biomc'
endif

use_python_define = '-DUSE_PYTHON_MODULE' # Define flag to pass to compiler if  use_dynamic_module is true

## Get buildtype specific options

use_kokkos_tools = get_option('use_kokkos_tools')
CMAKE_TYPE = 'Release'
if build_type == 'debug'
    # NOTHING
    CMAKE_TYPE = 'Debug'
    add_global_arguments('-DDEBUG', language: 'cpp')
    add_global_arguments('-DFIX_SEED', language: 'cpp')
elif build_type == 'debugoptimized'
    add_global_arguments('-DDEBUG', language: 'cpp')
    add_global_arguments('-DKOKKOS_ENABLE_DEBUG', language: 'cpp')
    if compiler.has_argument('-ffast-math')
        arch_specific_args = [arch_specific_args, '-ffast-math', '-Wno-nan-infinity-disabled']
    endif
else
    add_global_arguments('-DNDEBUG', language: 'cpp')

    if compiler.has_argument('-ffast-math')
        arch_specific_args = [arch_specific_args, '-ffast-math', '-Wno-nan-infinity-disabled']
    endif

endif

cpp_args = [arch_specific_args]
if (use_kokkos_tools)
    cpp_args = [cpp_args, '-DENABLE_KOKKOS_PROFILING']
endif

mpi = dependency('mpi', language: 'cpp', required: true, include_type: 'system')
thread_dep = dependency('threads', required: true, include_type: 'system')
omp = dependency('openmp', include_type: 'system')
hwloc = dependency('hwloc', include_type: 'system')
dl_dep = dependency('dl', required: true, include_type: 'system')

eigen_dep = dependency('eigen3', required: true, version: '>=5.0.1')
blas = dependency('blas', required: false) # Eigen can work without BLAS
cpp_args = [cpp_args, '-DEIGEN_NO_CUDA=1']

# add_global_arguments('-DEIGEN_DONT_VECTORIZE=1', language: 'cpp')
# if blas.found()
#   add_global_arguments('-DEIGEN_USE_BLAS', language: 'cpp')
# endif
#add_global_arguments('-DEIGEN_DONT_PARALLELIZE',language:'cpp')

cmake = import('cmake', required: true)
if use_system_kokkos
    __cmake_kokkos_dep__ = dependency(
        'Kokkos',
        method: 'cmake',
        required: true,
        components: ['kokkos'],
        modules: ['Kokkos::kokkos'],
    )

    kokkos_compile_arg = ['-DKOKKOS_DEPENDENCE', '-Wno-unknown-pragmas']

    if use_cuda

        kokkos_compile_arg = [
            kokkos_compile_arg,
            '-extended-lambda',
            '-Wext-lambda-captures-this',
            '--expt-relaxed-constexpr',
        ]

    endif

    kokkos = declare_dependency(
        dependencies: [__cmake_kokkos_dep__, omp, hwloc, thread_dep, dl_dep],
        compile_args: kokkos_compile_arg,
    )
else

    kokkos_option = cmake.subproject_options()
    kokkos_option.add_cmake_defines({'CMAKE_CXX_STANDARD': '20'})
    kokkos_option.add_cmake_defines({'CMAKE_BUILD_TYPE': CMAKE_TYPE})
    kokkos_option.add_cmake_defines({'Kokkos_ENABLE_OPENMP': 'ON'})
    kokkos_option.add_cmake_defines({'Kokkos_ENABLE_HWLOC': 'ON'})
    kokkos_option.add_cmake_defines({'CMAKE_POSITION_INDEPENDENT_CODE': 'ON'})
    kokkos_option.add_cmake_defines({'Kokkos_ENABLE_DEPRECATED_CODE_4': 'OFF'})
    kokkos_option.add_cmake_defines({'Kokkos_ENABLE_IMPL_VIEW_LEGACY': 'OFF'})
    kokkos_option.set_install(false)
    sub_kokkos = cmake.subproject('kokkos', required: true, options: kokkos_option)

    kokkos = declare_dependency(
        dependencies: [
            sub_kokkos.dependency('kokkos'),
            sub_kokkos.dependency('kokkoscore'),
            sub_kokkos.dependency('kokkossimd'),
            sub_kokkos.dependency('kokkoscontainers'),
            sub_kokkos.dependency('kokkosalgorithms'),
            omp,
            hwloc,
            thread_dep,
            dl_dep,
        ],
    )
endif

kokkos_sampling = dependency('kokkos_sampling-1.0', required: true)
# magic_enum= dependency('magic_enum', required: true)

if build_udf
    dynlib_dep = dependency('dynlib-1.0', required: true)
endif
highfive_found = false

ceral_found = false
if cmake.found()

    highfive_cmake_options = cmake.subproject_options()
    highfive_cmake_options.add_cmake_defines({'CMAKE_CXX_STANDARD': '20'})
    highfive_cmake_options.add_cmake_defines({'CMAKE_BUILD_TYPE': CMAKE_TYPE})
    highfive_cmake_options.add_cmake_defines({'HIGHFIVE_USE_BOOST': false})
    highfive_cmake_options.add_cmake_defines({'HIGHFIVE_EXAMPLES': false})
    highfive_cmake_options.add_cmake_defines({'HIGHFIVE_UNIT_TESTS': false})
    highfive_cmake_options.add_cmake_defines({'HIGHFIVE_USE_BOOST': false})
    ## MAYBE USE ?? IDK
    highfive_cmake_options.add_cmake_defines({'HDF5_PREFER_PARALLEL': false})
    highfive_cmake_options.add_cmake_defines({'HDF5_IS_PARALLEL': false})
    highfive_cmake_options.set_install(false)
    highfive_proj = cmake.subproject('highfive', options: highfive_cmake_options, required: false)

    if highfive_proj.found()
        hdf5 = dependency('hdf5', language: 'cpp', required: false)
        if hdf5.found()
            highfive_found = true
            add_project_arguments('-DUSE_HIGHFIVE', language: 'cpp')
            highfive_dependency = [
                highfive_proj.dependency('HighFive'),
                hdf5,
                highfive_proj.dependency('HighFiveInclude'),
            ]
        endif
    endif
    if use_cereal
        cereal_cmake_options = cmake.subproject_options()
        cereal_cmake_options.add_cmake_defines({'SKIP_PERFORMANCE_COMPARISON': true})
        cereal_cmake_options.add_cmake_defines({'CMAKE_BUILD_TYPE': CMAKE_TYPE})
        cereal_cmake_options.add_cmake_defines({'BUILD_DOC': false})
        cereal_cmake_options.add_cmake_defines({'BUILD_SANDBOX': false})
        cereal_cmake_options.add_cmake_defines({'BUILD_SANDBOX': false})
        cereal_cmake_options.add_cmake_defines({'BUILD_TESTS': false})
        cereal_cmake_options.set_install(false)
        cereal_proj = cmake.subproject('cereal', options: cereal_cmake_options, required: false)
        if cereal_proj.found()
            #and not use_cuda
            ceral_found = true
            if compiler.has_argument('-Wno-missing-template-arg-list-after-template-kw')
                cpp_args = [
                    cpp_args,
                    '-Wno-missing-template-arg-list-after-template-kw',
                ] #FIXME needed with clang++-19
            endif
            add_project_arguments('-DUSE_CEAREAL', language: 'cpp')
            cereal_dependency = cereal_proj.dependency('cereal')
        endif
    endif

endif

cma_dep = dependency('rcmtool', required: true, version: '>=0.1.0', static: true)

add_project_arguments(cpp_args, language: 'cpp')
subdir('apps')

_r = run_command('sh', '-c', f'ln -s @outdir@/compile_commands.json .', check: false)

subdir('tools')
if not python_install
    meson.add_install_script('devutils/install_script.sh', result_dir)
endif


message('select model(s)',get_option('model_list_name'))
if (use_kokkos_tools)
    message('Kokkos tool enable')
endif
