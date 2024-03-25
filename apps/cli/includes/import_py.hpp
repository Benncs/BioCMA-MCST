#ifndef __PY_MODULES_HPP__
#define __PY_MODULES_HPP__

#ifdef BIO_DYNAMIC_MODULE

#  include <simulation/models/types.hpp>
#  include <string>
#  include <pybind11/embed.h>

KModel get_python_module(const std::string &module_path);
pybind11::scoped_interpreter init_dynamic_module();

#else
#  error "get_python_module is onyl available with dynamic module activated"
#endif

#endif //__PY_MODULES_HPP__