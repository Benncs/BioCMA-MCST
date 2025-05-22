

test_probes = executable(
  'test_probes',
  'test_probes.cpp',
  dependencies: [simulation_lib_dependency],
  include_directories: [private_simulation_includes],
)


test_feed= executable(
  'test_feed',
  'test_feed.cpp',
  dependencies: [simulation_lib_dependency],
  include_directories: [private_simulation_includes],
)

# test_log = executable(
#   'test_log',
#   'test_log.cpp',
#   dependencies: [simulation_lib_dependency],
#   include_directories: [],
# )
# test('test_log', test_log)

test('test_probes', test_probes)
test('test_feed', test_feed)