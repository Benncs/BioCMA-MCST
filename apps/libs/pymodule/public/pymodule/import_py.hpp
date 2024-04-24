#ifndef __PY_MODULES_HPP__
#define __PY_MODULES_HPP__

#include <models/types.hpp>
#include <pybind11/embed.h>
#include <string>

KModel get_python_module(const std::string &module_path);
pybind11::scoped_interpreter init_python_interpreter();

#endif //__PY_MODULES_HPP__