t1 = executable('test_main', 'test_main.cpp',dependencies:[app_test_dependency])
t2 = executable('test_runtime_init', 'test_rt_init.cpp',dependencies:[app_test_dependency])
# test_scalar_factory = executable('test_scalar_factory', 'test_scalar_factory.cpp',dependencies:[app_test_dependency])

test('test main', t1)
test('test runtime_init', t2)


# if get_option('use_dynamic_module')
#   t3 = executable('test_test_pymodule', 'test_pymodule.cpp',dependencies:[app_test_dependency_dyn])
#   test('test test_pymodule',t3)
# endif 