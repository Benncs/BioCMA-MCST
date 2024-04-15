t1 = executable('test_probabilites', 'test_transport.cpp',dependencies:[simu_dependency],include_directories:['../src/'])
test_flow_map = executable('test_flowmap', 'test_flowmap.cpp',dependencies:[simu_dependency])


test('test_flowmap', test_flow_map)
test('test test_probabilites', t1)

