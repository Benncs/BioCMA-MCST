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

pybind11::scoped_interpreter init_python_interpreter()
{
  std::cout << "PYTHON INTERPRETER INITIALISATION" << std::endl;
  pybind11::scoped_interpreter interpret;

  return interpret;
}

model_properties_t convert_to_variant(const pybind11::handle& obj)
{
  if (py::isinstance<py::int_>(obj))
  {
    return obj.cast<int>();
  }
  else if (py::isinstance<py::float_>(obj))
  {
    // Check if the float is narrow enough to fit into a float, otherwise use
    // double
    double d = obj.cast<double>();
    if (d == static_cast<float>(d))
    {
      return static_cast<float>(d);
    }
    else
    {
      return d;
    }
  }
  else if (py::isinstance<py::str>(obj))
  {
    return obj.cast<std::string>();
  }
  else
  {
    throw std::invalid_argument(
        "Unsupported type for conversion to numeric_model_properties_t");
  }
}

model_properties_detail_t convert_dict_to_map(const py::dict &kwargs)
{
  std::unordered_map<std::string, model_properties_t> cpp_map;

  for (auto item : kwargs)
  {
    // Assuming the keys are strings
    std::string key = py::str(item.first);
    model_properties_t value = convert_to_variant(item.second);

    cpp_map[key] = value;
  }

  return cpp_map;
}

void init_module_function(KModel &model,
                          py::function &&init,
                          py::function &&update,
                          py::function &&contrib,
                          py::function &&division,
                          py::function &&properties,
                          py::function &&__debug)
{
  model.init_kernel = [init = std::move(init)](MC::Particles &p)
  {
    p.data = std::make_shared<OpaquePointer>();
    {
      init(std::ref(p));
      return;
    }
  };

  model.update_kernel =
      [update = std::move(update)](double dt, MC::Particles &p, auto &&c)
  {
    update(dt, std::ref(p), std::vector(c.begin(), c.end()));
    return;
  };

  model.contribution_kernel = [contrib = std::move(contrib)](auto &&p, auto &&o)
  {
    contrib(std::ref(p), std::ref(o));

    return;
  };
  model.division_kernel = [division = std::move(division)](auto &&p)
  {
    MC::Particles child(p);
    auto s = child.status;
    division(std::ref(p), std::ref(child));
    auto s2 = child.status;
    // assert(s != s2);
    return child;
  };

  model.get_properties =
      [properties = std::move(properties)](auto &&p)
  {
    auto dict = properties(std::cref(p));
    return convert_dict_to_map(dict);
  };
}

KModel get_python_module(const std::string &module_path)
{
  KModel model;
  try
  {
    py::module model_setup = py::module_::import(module_path.c_str());
    py::function init = model_setup.attr("init_kernel");
    py::function update = model_setup.attr("update_kernel");
    py::function contrib = model_setup.attr("contribution_kernel");
    py::function division = model_setup.attr("division_kernel");
    py::function properties = model_setup.attr("get_properties");
    py::function __debug = model_setup.attr("__debug");

    init_module_function(model,
                         std::move(init),
                         std::move(update),
                         std::move(contrib),
                         std::move(division),
                         std::move(properties),
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
