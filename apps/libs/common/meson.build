
src_dir = './src'

common_sources = run_command('sh', '-c', f'find @src_dir@ -type f -name *.cpp ').stdout().splitlines()

public_common = include_directories('public')

incdirs = [public_common, 'includes',auto_generate_include]

public_includes = [public_common,auto_generate_include]

common_lib = static_library(
   f'@project_name@_common',
   common_sources,
   include_directories: incdirs,
   dependencies: [kokkos],
)

common_dependecy = declare_dependency(
   dependencies:[kokkos],
   include_directories: public_includes,
   link_with: common_lib,
   link_args: [flto_flag]
)


if build_test
subdir('tests')
endif
