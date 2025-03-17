

if ceral_found  
test_serde_kview = executable('test_serde_kview', 'test_serde_kview.cpp',dependencies:[mc_dependency,cereal_dependency])
test_serde_container = executable('test_serde_container', 'test_serde_container.cpp',dependencies:[mc_dependency,cereal_dependency])
test('test_serde_kview',test_serde_kview)
test('test_serde_container',test_serde_container)
endif
# test_rng = executable('test_rng', 'test_rng.cpp',dependencies:[mc_dependency])
# test('test_rng',test_rng)



dl_dep = compiler.find_library('dl', required : false)
test_list = executable('test_list', 'test_list.cpp',dependencies:[mc_dependency,dl_dep])
test('test_list',test_list)

arg = []
if ci_execution
arg=['-DCI_TEST']
endif
test_rng_2 = executable('test_rng_2', 'test_rng_2.cpp',dependencies:[mc_dependency],cpp_args:arg)
test('test_rng_2',test_rng_2)

if not ci_execution
t2 = executable('test_rng_long', 'test_rng_long.cpp',dependencies:[mc_dependency])
endif 