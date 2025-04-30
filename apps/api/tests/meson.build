cma_path = meson.current_source_dir() + '/data/0d/'


test_api = executable('test_api', 'test_api.cpp',dependencies:[lib_api_shared_dep])
test_api_bad_alloc = executable('test_api_bad_alloc', 'test_api_throw.cpp',dependencies:[lib_api_shared_dep])
if build_c_api
test_api_raw = executable('test_api_raw', 'test_api_raw.cpp',dependencies:[lib_api_c_shared_dep])
test('api_api_raw',test_api_raw,args:[cma_path])
endif 
test('api',test_api,args:[cma_path])
test('api_bad_alloc',test_api_bad_alloc,args:[cma_path])
