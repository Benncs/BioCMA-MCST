

src_dir = './src'

app_sources = run_command('sh','-c', f'find @src_dir@ -type f -name *.cpp ').stdout().splitlines()


lib_sources = []
foreach item : app_sources
  if item not in [f'@src_dir@/main.cpp',f'@src_dir@/import_py.cpp']
    lib_sources += item
  endif
endforeach

message(lib_sources)

mpi=dependency('mpi', language: 'cpp') 


dep = [mpi,biomodel_dependency, omp,cma_dep,mc_dependency,simu_dependency,common_dependecy,thread_dep]
name = f'@project_name@_cli_app'


# Define static library and its dependency
lib_app = static_library(name,
                         lib_sources,
                         include_directories: ['./includes'],
                         dependencies: dep,
                         cpp_args: ['-UBIO_DYNAMIC_MODULE'])

app_test_dependency = declare_dependency(
    include_directories: ['./includes'],
    link_with: lib_app,
    dependencies: dep
)

# Define executable
executable(name,
                     f'@src_dir@/main.cpp',
                     include_directories: [ './includes'],
                     dependencies: [app_test_dependency],
                     cpp_args: ['-UBIO_DYNAMIC_MODULE'],cpp_args:target_custom_options)

# Conditional dynamic module build
if get_option('use_dynamic_module')

  name = name + '_dyn'

  lib_app_dyn = static_library(name,
                               lib_sources,
                               include_directories: [ './includes'],
                               dependencies: dep)

  app_test_dependency_dyn = declare_dependency(
      include_directories: ['./includes'],
      link_with: lib_app_dyn,
      dependencies: [dep,py_dynamic_module],
      sources: lib_sources + [f'@src_dir@/import_py.cpp']
  )

  # Assuming py_dynamic_module is a valid dependency
  executable(name,
                           f'@src_dir@/main.cpp',
                           include_directories: ['./includes'],
                           dependencies: [app_test_dependency_dyn, py_dynamic_module])
endif


subdir('tests')

