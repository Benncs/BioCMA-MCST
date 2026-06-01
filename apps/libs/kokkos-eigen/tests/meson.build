lib_app = library(
    f'@project_name@',
    include_directories: [],
    dependencies: [kokkos_eigen, kokkos, eigen_dep],
    cpp_args: [],
)

dep_app = declare_dependency(
    include_directories: [],
    link_with: lib_app,
    dependencies: [kokkos_eigen, kokkos, eigen_dep],
    compile_args: [],
)

executable(
    f'@project_name@',
    f'main.cpp',
    dependencies: [dep_app],
    link_args: [flto_flag],
)

flags_test = '-Wall -Wextra '

test_basic = executable(
    'kokkos_eigen_test_basic',
    'test_basic.cpp',
    dependencies: [dep_app],
    cpp_args: [flags_test],
)

test_sync = executable(
    'kokkos_eigen_test_sync',
    'test_sync.cpp',
    dependencies: [dep_app],
    cpp_args: [flags_test],
)
test('kokkos_eigen_test_basic', test_basic)
test('kokkos_eigen_test_sync', test_sync)