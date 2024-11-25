#ifndef __IMPORT_PY_HPP__
#define __IMPORT_PY_HPP__
#include <mc/particles/particle_model.hpp>
#include <memory>

#define EXPORT __attribute__((visibility("default")))

namespace pybind11
{
  class scoped_interpreter;
} // namespace pybind11


using python_interpreter_t =
    std::unique_ptr<pybind11::scoped_interpreter, void (*)(pybind11::scoped_interpreter*)>;

python_interpreter_t init_python_interpreter();

#endif
