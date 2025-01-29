py = import('python').find_installation(pure: false)
pybind11_dep = dependency('pybind11', required: true)
pydep = [py.dependency(embed: true),pybind11_dep]



executable(
    'model_fuzzing',
    ['main.cpp','embed.cpp'],
    dependencies: [biomodel_dependency,kokkos,core_shared_dependency,pydep],
    include_directories: [public_include_dir],
)
