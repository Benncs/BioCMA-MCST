
src_dir = './src'

mc_sources = run_command('sh', '-c', f'find @src_dir@ -type f -name *.cpp ').stdout().splitlines()

public_mc = include_directories('public')

public_include_dir = [public_mc, public_common,models_list_def,auto_generate_include]

arg = []
if use_dynamic_module
public_include_dir=[public_include_dir,py_wrap_inc]
arg=[use_python_define]
endif 

mc_lib = static_library(
   f'@project_name@_mc',
   mc_sources,
   include_directories: [public_include_dir, 'includes'],
   dependencies: [kokkos],
   # pic: true,
   # install: true,
   install: false,
   install_rpath: global_install_dir,
   link_args: [flto_flag],
    cpp_args: [arg],
)

mc_dependency = declare_dependency(
   include_directories: public_include_dir,
   dependencies:[kokkos],
   link_with: mc_lib,
   compile_args: [arg],
)

subdir('./tests')

# install_subdir(
#    'public',
#    install_dir: includes_install_dir,
# )

