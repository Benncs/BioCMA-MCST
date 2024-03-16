

src_dir = './src'

app_sources = run_command('sh','-c', f'find @src_dir@ -type f -name *.cpp ').stdout().splitlines()


lib_sources = []
foreach item : app_sources
  if item not in [f'@src_dir@/main.cpp']
    lib_sources += item
  endif
endforeach

message(lib_sources)

cma_lib_dir = '/home/benjamin/Documenti/code/cpp/cma_read/builddir'
cma_public_dir = '/home/benjamin/Documenti/code/cpp/cma_read/public'

cma_dep = cpp.find_library('cma_read_static', dirs:cma_lib_dir)




mpi=dependency('mpi', language: 'cpp') 
omp = dependency('openmp',language:'cpp')



lib_app = static_library(f'@project_name@_cli_app',lib_sources,include_directories : [cma_public_dir,'./includes'],
                            dependencies : [mpi, omp,cma_dep,mc_dependency,simu_dependency,common_dependecy])


app_test_dependency  = declare_dependency(
   include_directories : [cma_public_dir,'./includes'], 
   link_with : lib_app,
   dependencies:[mpi, omp,cma_dep,mc_dependency,simu_dependency,common_dependecy],
   sources : lib_sources)

app_cli = executable(f'@project_name@_cli_app',
                            f'@src_dir@/main.cpp',
                            include_directories : [cma_public_dir,'./includes'],
                            dependencies : [app_test_dependency])
                            



subdir('tests')

