
public_includes = include_directories(f'public/')
kokkos_eigen = declare_dependency(
    include_directories: public_includes,
    compile_args:['-DEIGEN_NO_CUDA=1']
)
if build_test
subdir('tests')
endif
