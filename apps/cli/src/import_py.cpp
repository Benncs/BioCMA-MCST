#ifdef BIO_DYNAMIC_MODULE
#  include "pymodule/opaque_type.hpp"

#  include <import_py.hpp>
#  include <pybind11/embed.h>
#  include <pybind11/gil.h>
#  include <pybind11/numpy.h>
#  include <pybind11/pytypes.h>
#  include <pybind11/stl.h>
#  include <pymodule/opaque_type.hpp>
#  include "mc/particles/mcparticles.hpp"
#  include <omp.h>
#  include <simulation/models/types.hpp>
#  include <stdexcept>
#  include <memory>

namespace py = pybind11;
static PyGILState_STATE gstate;


pybind11::scoped_interpreter init_dynamic_module()
{
  pybind11::scoped_interpreter interpret;
  return interpret;
}



KModel get_python_module(const std::string &module_path)
{
  KModel model;
  try
  {

    py::module model_setup = py::module::import(module_path.c_str());

    py::function init = model_setup.attr("init_kernel");
    py::function update = model_setup.attr("update_kernel");
    model.init_kernel = [init](MC::Particles &p)
    {
      p.data = std::make_shared<OpaquePointer>();
      {

        try
        {
#  pragma omp critical
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
        [update = std::move(update)](
            double dt, MC::Particles &p, std::span<const double> c)
    {
      {
        try
        {
#  pragma omp critical
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

    model.contribution_kernel = [](auto &&p, auto &&o) { return; };
    model.division_kernel = [](auto &&p) { return; };
    return model;
  }

  catch (...)
  {
    throw std::runtime_error("Error loading python module");
  }
}

#endif