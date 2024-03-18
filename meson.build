project('BioCMA-MCST','cpp','c',version:'0.0.1',default_options:['warning_level=3','cpp_std=gnu++20','b_colorout=always'])



project_name = meson.project_name()
build_type = get_option('buildtype')

default_options = {
    'debug': ['b_sanitize=thread'],
    'release': ['b_lto=true'],
}

# Get default options based on the build type
default_opts = default_options.get(build_type, [])

# Add default options
foreach opt: default_opts
    add_project_arguments(opt, language: 'c')
endforeach


thread_dep = dependency('threads')
cpp = meson.get_compiler('cpp')

eigen_includes = include_directories('./subprojects/eigen')

subdir('apps')

run_command('sh','-c','ln -s builddir/compile_commands.json .')


root = meson.project_source_root()
run_target('run_shared',command:app_cli)
run_target('bench_shared',command:['/bin/bash','-c', f'@root@/devutils/bench.sh '+ app_cli.full_path()])
run_target('bench_distributed',command:['/bin/bash','-c', f'@root@/devutils/bench.sh mpiexec '+ app_cli.full_path()])


run_target('run_distributed',command:['/bin/bash','-c','mpiexec '+ app_cli.full_path()])

benchmark('bench_shared', find_program('./devutils/bench.sh'), args: [app_cli.full_path()],timeout:180)
