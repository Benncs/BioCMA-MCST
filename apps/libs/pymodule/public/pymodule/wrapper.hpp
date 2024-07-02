#ifndef __PY_WRAPPER_HPP__
#define __PY_WRAPPER_HPP__

#include <pybind11/embed.h>
#include <memory>

class PythonWrapper
{
public:
explicit PythonWrapper(const std::string &module_path);
private:
    pybind11::scoped_interpreter m_interpreter;

    std::unique_ptr<pybind11::gil_scoped_release> mp_gil_release;
};



#endif 