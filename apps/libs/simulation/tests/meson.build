
test_transport = executable(
  'test_probabilites',
  'test_transport.cpp',
  dependencies: [simulation_lib_dependency],
  include_directories: ['../src/', private_simulation_includes],
)
test_flow_map = executable(
  'test_flowmap',
  'test_flowmap.cpp',
  dependencies: [simulation_lib_dependency],
  include_directories: [private_simulation_includes],
)
test_transition_matrix = executable(
  'test_transition_matrix',
  'test_transition_matrix.cpp',
  dependencies: [simulation_lib_dependency],
  include_directories: [private_simulation_includes],
)


test('test_flowmap', test_flow_map)
test('test test_probabilites', test_transport)
test('test transition matrix', test_transition_matrix)
