dep = [mc_dependency, dynlib_dep]
foreach name : ['minimal', 'division_time']
  library(
    f'udf_model_@name@',
    f'@name@.cpp',
    include_directories: [include_directories('./'), public_model, ext_include],
    dependencies: dep,
    link_args: [flto_flag],
    pic: true,
    gnu_symbol_visibility: 'hidden',
    cpp_args: ['-DDECLARE_EXPORT_UDF'],
  )
endforeach