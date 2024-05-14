
src_dir = './src'

model_sources = run_command('sh', '-c', f'find @src_dir@ -type f -name *.cpp ').stdout().splitlines()

public_model = include_directories('public')


public_include_dir = [public_model]

dep = [mc_dependency, blas, common_dependecy, eigen_dep]

model_lib = static_library(
   f'@project_name@_bmodels',
   model_sources,
   include_directories: [public_include_dir, 'includes'],
   dependencies: dep,
   link_args: ['-flto'],
)

biomodel_dependency = declare_dependency(
   include_directories: [public_include_dir],
   link_with: model_lib,
   dependencies: dep,
)

install_subdir(
   'public',
   install_dir: includes_install_dir,
)

subdir('tests')
