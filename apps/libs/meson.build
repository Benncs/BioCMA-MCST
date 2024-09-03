
prog_python = find_program('python3')

# custom target to generate model_loader.cpp and model_list.hpp
generated_model_selection = custom_target(
    'generated_model_selection.[ch]',
    output: ['wrap_init_model_selector.cpp', 'wrap_init_model_selector.hpp', 'variant_model.hpp'],
    input: f'@root@/devutils/auto_generate/generate_model_loader.py',
    command: [
        prog_python,
        '@INPUT@',
        f'@root@/apps/libs/biomodels',
        f'@root@/devutils/auto_generate/model_loader.cxx',
        '@OUTPUT0@',
        f'@root@/devutils/auto_generate/model_loader.hxx',
        '@OUTPUT1@',
        f'@root@/devutils/auto_generate/variant_model.hxx',
        '@OUTPUT2@',
    ],
)

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
message(cpp_args)

#Static library that includes the generated files
model_selection = library(
    'model_selection',
    [generated_model_selection],
    include_directories: 'mc/public/',
    dependencies: dep,
    cpp_args: [cpp_args],
)

model_selection_dep = declare_dependency(
    sources: [generated_model_selection[1]],
    link_with: [model_selection],
    dependencies: dep,
    include_directories: 'mc/public/',
    compile_args: [cpp_args],
)