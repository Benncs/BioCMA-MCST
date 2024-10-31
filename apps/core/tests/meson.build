test_scalar_factory = executable('test_scalar_factory', 'test_scalar_factory.cpp',dependencies:[core_shared_dependency])


test_data_path = meson.current_source_dir() + '/test_scalar_read.h5'


test('test_scalar_factory',test_scalar_factory,args:[test_data_path])