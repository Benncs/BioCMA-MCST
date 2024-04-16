project('BioCMA-MCST','cpp',version:'0.0.2',default_options:['warning_level=3','cpp_std=gnu++20','b_colorout=always'])


root = meson.source_root()
project_name = meson.project_name()
build_type = get_option('buildtype')
cpp = meson.get_compiler('cpp')

default_options = {
    # 'debug': ['b_sanitize=thread'],
    # 'release': ['b_lto=true'],
}

# Get default options based on the build type
default_opts = default_options.get(build_type, [])

if build_type == 'debug'
    add_global_arguments('-DDEBUG', language : 'cpp')
endif

foreach opt: default_opts
    add_project_arguments(opt, language: 'cpp')
endforeach


thread_dep = dependency('threads')
omp = dependency('openmp',language:'cpp',required:true)


cpp_args = ['-march=native','-fno-math-errno']

add_project_arguments(cpp_args, language : 'cpp')

blas = dependency('blas',required:false)
eigen_dep = dependency('eigen3')
if blas.found()
add_global_arguments('-DEIGEN_USE_BLAS', language : 'cpp')
endif

use_dynamic_module = get_option('use_dynamic_module')



if use_dynamic_module
    message('Dynamic module feature activated')
    pybind11_dep = dependency('pybind11', required : true)
    if pybind11_dep.found()
    add_global_arguments('-DBIO_DYNAMIC_MODULE=ON', language : 'cpp')
        python_embded = declare_dependency(
            include_directories : include_directories('/usr/include/python3.10', pybind11_dep.get_pkgconfig_variable('includedir')),
            dependencies : [pybind11_dep],
            link_args : ['-L/usr/lib', '-lpython3.10']
        )
    endif
endif 


cma_reader = subproject('birem',default_options:['warning_level=3','cpp_std=gnu++20','b_colorout=always','python_wrap=false'])

cma_dep = cma_reader.get_variable('cma_read_shared_lib_dependency') 

subdir('apps')

run_command('sh','-c','ln -s builddir/compile_commands.json .')


root = meson.project_source_root()
# run_target('run_shared',command:app_cli)
# run_target('bench_shared',command:['/bin/bash','-c', f'@root@/devutils/bench.sh '+ app_cli.full_path()])
# run_target('bench_distributed',command:['/bin/bash','-c', f'@root@/devutils/bench.sh mpiexec '+ app_cli.full_path()])


# run_target('run_distributed',command:['/bin/bash','-c','mpiexec '+ app_cli.full_path()])

# benchmark('bench_shared', find_program('./devutils/bench.sh'), args: [app_cli.full_path()],timeout:180)

