
models_list_def = include_directories('./biomodels/public')


subdir('common')
subdir('mpi_w')
subdir('mc')
subdir('biomodels')

if use_dynamic_module
    subdir('pymodule')
endif

subdir('simulation')


dep = [biomodel_dependency]
cpp_args = []
if use_dynamic_module
    dep = [biomodel_dependency, py_dynamic_module]
    cpp_args = [use_python_define]
endif

#Static library that includes the generated files
model_selection = static_library(
    'model_selection',
    [generated_model_selection],
    include_directories: public_mc,
    dependencies: dep,
    cpp_args: [cpp_args],
)

model_selection_dep = declare_dependency(
    sources: [generated_model_selection[1]],
    link_with: [model_selection],
    dependencies: dep,
    include_directories: public_mc,
    compile_args: [cpp_args],
)