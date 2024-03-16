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



cpp = meson.get_compiler('cpp')
subdir('apps/libs')
subdir('apps/cli')

run_command('sh','-c','ln -s builddir/compile_commands.json .')