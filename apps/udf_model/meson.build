

dep = [mc_dependency]

library(
    'udf_model',
   'impl.cpp',
   include_directories: [include_directories('./'),shared_load_include,public_model,ext_include],
   dependencies: dep,
   link_args: [flto_flag],
   pic:true,
   gnu_symbol_visibility: 'hidden',
   cpp_args:['-DDECLARE_EXPORT_UDF']
)