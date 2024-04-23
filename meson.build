project(
  'BioCMA-MCST',
  'cpp',
  version: '0.0.2',
  default_options: [
    'warning_level=3',
    'cpp_std=gnu++20',
    'b_colorout=always',
    'b_coverage=true',
  ],
)


root = meson.project_source_root()
project_name = meson.project_name()
build_type = get_option('buildtype')
cpp = meson.get_compiler('cpp')


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



use_dynamic_module = get_option('use_dynamic_module')
use_python_define = '-DUSE_PYTHON_MODULE'
cpp_args = ['-march=native', '-fno-math-errno']

mpi = dependency('mpi', language: 'cpp', required: true)
thread_dep = dependency('threads', required: true)
omp = dependency('openmp', language: 'cpp', required: true)
eigen_dep = dependency('eigen3', required: true)

blas = dependency('blas', required: false)
if blas.found()
  add_global_arguments('-DEIGEN_USE_BLAS', language: 'cpp')
endif


add_project_arguments(cpp_args, language: 'cpp')



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
  ],
)

cma_dep = cma_reader.get_variable('cma_read_shared_lib_dependency')

subdir('apps')

run_command('sh', '-c', 'ln -s builddir/compile_commands.json .', check: false)


# root = meson.project_source_root()

# clangfiy = find_program('clang-format', required: false)

# all_sources = run_command('sh','-c', f'find apps/ -type f -name *.cpp ').stdout().splitlines()


# run_target('run_shared',command:app_cli)
# run_target('bench_shared',command:['/bin/bash','-c', f'@root@/devutils/bench.sh '+ app_cli.full_path()])
# run_target('bench_distributed',command:['/bin/bash','-c', f'@root@/devutils/bench.sh mpiexec '+ app_cli.full_path()])


# run_target('run_distributed',command:['/bin/bash','-c','mpiexec '+ app_cli.full_path()])

# benchmark('bench_shared', find_program('./devutils/bench.sh'), args: [app_cli.full_path()],timeout:180)

