
t1 = executable('test_probabilites', 'test_transport.cpp',dependencies:[simu_dependency,simulation_public_dep],include_directories:['../src/',private_simu_includes])
test_flow_map = executable('test_flowmap', 'test_flowmap.cpp',dependencies:[simu_dependency,simulation_public_dep],include_directories:[private_simu_includes])
test_transition_matrix= executable('test_transition_matrix', 'test_transition_matrix.cpp',dependencies:[simu_dependency,simulation_public_dep],include_directories:[private_simu_includes])


test('test_flowmap', test_flow_map)
test('test test_probabilites', t1)
test('test transition matrix', test_transition_matrix)
