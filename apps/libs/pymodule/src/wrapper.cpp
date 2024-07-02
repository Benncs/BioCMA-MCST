#include <pymodule/wrapper.hpp>

// class PythonWrapper {
// public:
//     PythonWrapper(const std::string &module_path) : m_interpreter() {
//         py::object obj = py::module::import("main").attr("PythonClass")();  
        
//         mp_gil_release = std::make_unique<py::gil_scoped_release>();
//     }
// private:
//     py::scoped_interpreter m_interpreter;

//     std::unique_ptr<py::gil_scoped_release> mp_gil_release;
// };

PythonWrapper::PythonWrapper(const std::string &module_path)
{
    pybind11::object obj = pybind11::module::import(module_path.c_str());
    mp_gil_release = std::make_unique<pybind11::gil_scoped_release>();
}