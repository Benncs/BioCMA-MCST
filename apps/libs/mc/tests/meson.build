test_container = executable('test_container', 'test_container.cpp',dependencies:[mc_dependency])
test('test_container',test_container)


if ceral_found  and not use_cuda
test_serde_kview = executable('test_serde_kview', 'test_serde_kview.cpp',dependencies:[mc_dependency,cereal_dependency])
test_serde_container = executable('test_serde_container', 'test_serde_container.cpp',dependencies:[mc_dependency,cereal_dependency])
test('test_serde_kview',test_serde_kview)
test('test_serde_container',test_serde_container)
endif

arg = []
if ci_execution
arg=['-DCI_TEST']
endif
test_rng_2 = executable('test_rng_2', 'test_rng_2.cpp',dependencies:[mc_dependency],cpp_args:arg)
test('test_rng_2',test_rng_2,timeout: -1)

# if not ci_execution
# t2 = executable('test_rng_long', 'test_rng_long.cpp',dependencies:[mc_dependency])
# test('test_rng_long',t2,timeout: -1)
# endif 

