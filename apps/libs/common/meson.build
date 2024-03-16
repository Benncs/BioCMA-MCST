
src_dir = './src'

common_sources = run_command('sh','-c', f'find @src_dir@ -type f -name *.cpp ').stdout().splitlines()

public_common  =  include_directories('public')

incdirs = [public_common,'includes']

common_lib = library(f'@project_name@_common', common_sources,include_directories :incdirs,pic:true)

common_dependecy = declare_dependency(
   include_directories : incdirs,
   link_with : common_lib,
   sources : common_sources)


