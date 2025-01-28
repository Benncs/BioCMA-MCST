#t2 = executable('test_rng_long', 'test_rng_long.cpp',dependencies:[mc_dependency])
#test_rng = executable('test_rng', 'test_rng.cpp',dependencies:[mc_dependency])

if ceral_found
#test_serde_particle = executable('test_serde_particle', 'test_serde_particle.cpp',dependencies:[mc_dependency,cereal_dependency])
#test_serde_list = executable('test_serde_list', 'test_serde_list.cpp',dependencies:[mc_dependency,cereal_dependency])
#test('test_serde_particle',test_serde_particle)
#test('test_serde_list',test_serde_list,cpp_args:['-ldl'])
endif
#test('test_rng',test_rng)

# test('test_serde_particle',test_serde_particle)


test_list = executable('test_list', 'test_list.cpp',dependencies:[mc_dependency])
test('test_list',test_list)