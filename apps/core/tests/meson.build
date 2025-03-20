# test_scalar_factory = executable('test_scalar_factory', 'test_scalar_factory.cpp',dependencies:[core_shared_dependency])

test_signal_handler = executable(
    'test_signal_handler',
    'test_signal_handler.cpp',
    dependencies: [core_shared_dependency],
    include_directories: private_core_includes,
)

test_postprocess = executable(
    'test_postprocess',
    'test_postprocess.cpp',
    dependencies: [core_shared_dependency],
    include_directories: private_core_includes,
)

test_load_balancing = executable(
    'test_load_balancing',
    'test_load_balancing.cpp',
    dependencies: [core_distributed_dependency],
    include_directories: private_core_includes,
)

test_data_path = meson.current_source_dir() + '/test_scalar_read.h5'
# test('core_scalar_factory',test_scalar_factory,args:[test_data_path]) # TODO FIX CMAREAD VIEW
test('core_signal_handler', test_signal_handler)
test('test_postprocess', test_postprocess)
test('test_load_balancing', test_load_balancing)