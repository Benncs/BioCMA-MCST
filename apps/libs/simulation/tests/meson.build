# test_transport = executable(
#   'test_probabilites',
#   'test_transport.cpp',
#   dependencies: [simulation_lib_dependency],
#   include_directories: ['../src/', private_simulation_includes],
# )
test_flow_map = executable(
  'test_flowmap',
  'test_flowmap.cpp',
  dependencies: [simulation_lib_dependency],
  include_directories: [private_simulation_includes],
)
# test_transition_matrix = executable(
#   'test_transition_matrix',
#   'test_transition_matrix.cpp',
#   dependencies: [simulation_lib_dependency],
#   include_directories: [private_simulation_includes],
# )

test_probes = executable(
  'test_probes',
  'test_probes.cpp',
  dependencies: [simulation_lib_dependency],
  include_directories: [private_simulation_includes],
)

test_log = executable(
  'test_log',
  'test_log.cpp',
  dependencies: [simulation_lib_dependency],
  include_directories: [],
)
# test('test_log', test_log)

test('test_probes', test_probes)
test('test_flowmap', test_flow_map)

# test('test test_probabilites', test_transport)
# test('test transition matrix', test_transition_matrix)