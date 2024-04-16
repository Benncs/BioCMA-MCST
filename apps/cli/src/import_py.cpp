#ifdef BIO_DYNAMIC_MODULE
#  include "pymodule/opaque_type.hpp"

#  include "mc/particles/mcparticles.hpp"
#  include <import_py.hpp>
#  include <memory>
#  include <models/types.hpp>
#  include <omp.h>
#  include <pybind11/eigen.h>
#  include <pybind11/embed.h>
#  include <pybind11/gil.h>
#  include <pybind11/numpy.h>
#  include <pybind11/pytypes.h>
#  include <pybind11/stl.h>
#  include <pymodule/opaque_type.hpp>
#  include <stdexcept>

namespace py = pybind11;
static PyGILState_STATE gstate;

pybind11::scoped_interpreter init_dynamic_module()
{
  pybind11::scoped_interpreter interpret;
  return interpret;
}

void init_module_function(KModel &model,
                          py::function init,
                          py::function update,
                          py::function contrib)
{
  model.init_kernel = [init](MC::Particles &p)
  {
    p.data = std::make_shared<OpaquePointer>();
    {

      try
      {
        // #  pragma omp critical
        {
          Py_BEGIN_ALLOW_THREADS;
          gstate = PyGILState_Ensure();
          init(std::ref(p));
          PyGILState_Release(gstate);
          Py_END_ALLOW_THREADS;
        }
      }
      catch (py::error_already_set &e)
      {

        e.discard_as_unraisable(__func__);
      }

      return;
    }
  };

  model.update_kernel =
      [update = std::move(update)](double dt, MC::Particles &p, auto &&c)
  {
    {
      try
      {
        // #  pragma omp critical
        {
          Py_BEGIN_ALLOW_THREADS;
          gstate = PyGILState_Ensure();
          update(dt, std::ref(p), std::vector(c.begin(), c.end()));
          PyGILState_Release(gstate);
          Py_END_ALLOW_THREADS;
        }
      }
      catch (py::error_already_set &e)
      {

        e.discard_as_unraisable(__func__);
      }

      return;
    }
  };

  model.contribution_kernel = [contrib](auto &&p, auto &&o)
  {
    {
      try
      {
        // #  pragma omp critical
        {
          Py_BEGIN_ALLOW_THREADS;
          gstate = PyGILState_Ensure();
          contrib(std::ref(p), std::ref(o));
          PyGILState_Release(gstate);
          Py_END_ALLOW_THREADS;
        }
      }
      catch (py::error_already_set &e)
      {

        e.discard_as_unraisable(__func__);
      }

      return;
    }
  };
  model.division_kernel = [](auto &&p) { return MC::Particles(); };
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
    init_module_function(model, init, update, contrib);
    return model;
  }

  catch (...)
  {
    throw std::runtime_error("Error loading python module");
  }
}

#endif