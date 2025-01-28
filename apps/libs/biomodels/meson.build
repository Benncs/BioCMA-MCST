src_dir = './src'

model_sources = run_command('sh', '-c', f'find @src_dir@ -type f -name model_*.cpp ').stdout().splitlines()

public_model = include_directories('public')

ext_include = include_directories('./ext')

public_include_dir = [public_model]

if not use_cuda

loader = static_library(
   'loader',
   f'ext/udf_includes.cpp',
   include_directories: [ ext_include,public_include_dir],
   dependencies: [mc_dependency, common_dependecy,dynlib_dep],
   cpp_args: ['-DDECLARE_EXPORT_UDF'],
)

udf_loader = declare_dependency(
   link_with: loader,
   dependencies: [mc_dependency, common_dependecy,dynlib_dep],
)
else
udf_loader = declare_dependency()
endif 

dep = [mc_dependency, common_dependecy, kokkos,udf_loader]

model_lib = static_library(
   f'@project_name@_bmodels',
   model_sources,
   include_directories: [public_include_dir, 'includes', ext_include],
   dependencies: dep,
   link_args: [flto_flag],
   # pic:true
)

biomodel_dependency = declare_dependency(
   include_directories: [public_include_dir, ext_include],
   link_with: model_lib,
   dependencies: dep,
)



# install_subdir(
#    'public',
#    install_dir: includes_install_dir,
# )

subdir('tests')
