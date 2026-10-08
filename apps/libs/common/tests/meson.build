t2 = executable('test_newton', 'test_newton.cpp', dependencies: [common_dependecy])
test_env_var = executable('test_env_var', 'test_env_var.cpp', dependencies: [common_dependecy])
test_format = executable('test_format', 'test_format.cpp', dependencies: [common_dependecy])
test_league_size = executable('test_league_size', 'test_league_size.cpp', dependencies: [common_dependecy])
test_maths = executable('test_maths', 'test_maths.cpp', dependencies: [common_dependecy])
test_traits = executable('test_traits', 'test_traits.cpp', dependencies: [common_dependecy])
test_results = executable('test_results', 'test_results.cpp', dependencies: [common_dependecy])



# test_kokkos_vector = executable('test_kokkos_vector', 'test_kokkos_vector.cpp',dependencies:[common_dependecy])
# test('test_kokkos_vector', test_kokkos_vector  )
test('test_newton', t2)
test('test_env_var', test_env_var)
test('test_format', test_format)
test('test_league_size', test_league_size)
test('test_maths', test_maths)
test('test_traits', test_traits)
test('test_results', test_results)
