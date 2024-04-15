
src_dir = './src'

common_sources = run_command('sh','-c', f'find @src_dir@ -type f -name *.cpp ').stdout().splitlines()

public_common  =  include_directories('public')

incdirs = [public_common,'includes']

public_includes = [public_common]

common_lib = library(f'@project_name@_common', common_sources,include_directories :incdirs,pic:true,dependencies:[omp])

common_dependecy = declare_dependency(sources:common_sources,
   include_directories : public_includes,
   link_with : common_lib)


subdir('tests')