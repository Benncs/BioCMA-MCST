src_dir = './src'

# model_sources = run_command('sh', '-c', f'find @src_dir@ -type f -name model_*.cpp ').stdout().splitlines()

model_sources = [f'@src_dir@/two_meta_div.cpp']

public_model = include_directories('public')

# ext_include = include_directories('./ext')
ext_include = include_directories('ext')
public_include_dir = [public_model]

if not use_cuda
model_sources = [model_sources,f'@src_dir@/udfmodel_user.cpp']
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
     compile_args:['-DDECLARE_EXPORT_UDF']
)

model_lib = static_library(
   f'@project_name@_bmodels',
   model_sources,
   include_directories: [public_include_dir],
   dependencies: [mc_dependency, common_dependecy],
   link_args: [flto_flag],
   # pic:true
)



else
udf_loader = declare_dependency()

endif 

dep = [mc_dependency, common_dependecy, kokkos,udf_loader]

if not use_cuda
biomodel_dependency = declare_dependency(
   link_with: model_lib,
   include_directories: [public_include_dir],
   dependencies: dep,
)
else
biomodel_dependency = declare_dependency(
   include_directories: [public_include_dir],
   dependencies: dep,
)
endif 


if build_test
subdir('tests')
endif
