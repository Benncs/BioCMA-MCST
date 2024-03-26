#include <import_py.hpp>
#include <iostream>



void test_get_python_module() {
    std::string module_path = "tests.test_pymodule";
    try {
        KModel model = get_python_module(module_path);
         assert(true);
    } catch (const std::runtime_error& e) {
       assert(false);
    }
}

int main()
{
  test_get_python_module();
}

