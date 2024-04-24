#include "pymodule/opaque_type.hpp"

#include "mc/particles/mcparticles.hpp"
#include <iostream>
#include <memory>
#include <models/types.hpp>
#include <omp.h>
#include <pybind11/eigen.h>
#include <pybind11/embed.h>
#include <pybind11/gil.h>
#include <pybind11/numpy.h>
#include <pybind11/pytypes.h>
#include <pybind11/stl.h>
#include <pymodule/import_py.hpp>
#include <sstream>
#include <stdexcept>

namespace py = pybind11;
static PyGILState_STATE gstate;

pybind11::scoped_interpreter init_python_interpreter()
{
  std::cout << "PYTHON INTERPRETER INITIALISATION" << std::endl;
  pybind11::scoped_interpreter interpret;
  return interpret;
}
// #define PYBIND11_NO_ASSERT_GIL_HELD_INCREF_DECREF

#define PY_LOCK_KRNL(__func_call__) __func_call__;

//                                      \
  // _Pragma("omp critical")                                                      \
  // {                                                                            \
  //   Py_BEGIN_ALLOW_THREADS gstate = PyGILState_Ensure();                       \
  //   __func_call__;                                                             \
  //   PyGILState_Release(gstate);                                                \
  //   Py_END_ALLOW_THREADS                                                       \
  // }

#define PY_TRY_LOCK_KRNL(__func_call__)                                        \
  {                                                                            \
    try                                                                        \
    {                                                                          \
      PY_LOCK_KRNL(__func_call__)                                              \
    }                                                                          \
    catch (py::error_already_set & e)                                          \
    {                                                                          \
      e.discard_as_unraisable(__func__);                                       \
    }                                                                          \
  }

void init_module_function(KModel &model,
                          py::function &&init,
                          py::function &&update,
                          py::function &&contrib,
                          py::function &&division,
                          py::function &&__debug)
{
  model.init_kernel = [init = std::move(init)](MC::Particles &p)
  {
    p.data = std::make_shared<OpaquePointer>();
    {
      PY_TRY_LOCK_KRNL(init(std::ref(p)));
      return;
    }
  };

  model.update_kernel =
      [update = std::move(update)](double dt, MC::Particles &p, auto &&c)
  {
    PY_TRY_LOCK_KRNL(update(dt, std::ref(p), std::vector(c.begin(), c.end())));
    return;
  };

  model.contribution_kernel = [contrib = std::move(contrib)](auto &&p, auto &&o)
  {
    PY_TRY_LOCK_KRNL(contrib(std::ref(p), std::ref(o)));

    return;
  };
  model.division_kernel = [division = std::move(division)](auto &&p)
  {
    MC::Particles child(p);
    auto s = child.status;
    PY_TRY_LOCK_KRNL(division(std::ref(p), std::ref(child)));
    auto s2 = child.status;
    assert(s != s2);
    return child;
  };
  // #ifdef DEBUG
  model.f_dbg = [__debug = std::move(__debug)](auto &&p)
  { PY_TRY_LOCK_KRNL(__debug(std::ref(p))); };
  // #endif
}

KModel get_python_module(const std::string &module_path)
{
  KModel model;
  try
  {
    py::module model_setup = py::module::import(module_path.c_str());
    py::function init = model_setup.attr("init_kernel");
    py::function update = model_setup.attr("update_kernel");
    py::function contrib = model_setup.attr("contribution_kernel");
    py::function division = model_setup.attr("division_kernel");
    py::function __debug = model_setup.attr("__debug");

    init_module_function(model,
                         std::move(init),
                         std::move(update),
                         std::move(contrib),
                         std::move(division),
                         std::move(__debug));
    return model;
  }
  catch (const std::exception &e)
  {
    std::stringstream msg;
    msg << "Error loading python module: " << std::endl;
    msg << e.what() << std::endl;
    throw std::runtime_error(msg.str());
  }
  catch (...)
  {
    throw std::runtime_error("Error loading python module");
  }
}
