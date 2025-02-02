#include "pybind11/attr.h"
#include "pybind11/buffer_info.h"
#include "pymodule/opaque_type.hpp"
#include <cstdlib>
#include <iostream>
#include <pybind11/embed.h>
#include <pybind11/gil.h>
#include <pybind11/numpy.h>
#include <pybind11/pytypes.h>
#include <pybind11/stl.h>
#include <pymodule/import_py.hpp>
#include <pymodule/model_python.hpp>

static bool interpreter_init = false;

static void _custom_interpreter_deleter(pybind11::scoped_interpreter* interpreter)
{
  if (interpreter != nullptr)
  {
    delete interpreter; // NOLINT
    interpreter_init = false;
    interpreter = nullptr;
    std::cout << "PYTHON INTERPRETER DESTROYED" << std::endl;
  }
}

EXPORT python_interpreter_t init_python_interpreter()
{
  if (!interpreter_init)
  {
    std::cout << "PYTHON INTERPRETER INITIALIZATION" << std::endl;
    interpreter_init = true;
    return {new pybind11::scoped_interpreter(), &_custom_interpreter_deleter};
  }
  return {nullptr, &_custom_interpreter_deleter};
}

// NOLINTBEGIN
#define PYTHON_CST_FWD(__value__) std::cref((__value__))
#define PYTHON_FWD(__value__) std::ref((__value__))
// NOLINTEND

namespace PythonWrap
{

  model_properties_detail_t convert_dict_to_map(const py::dict& kwargs)
  {
    model_properties_detail_t cpp_map;

    for (auto item : kwargs)
    {

      std::string key = py::str(item.first);
      cpp_map[key] = item.second.cast<double>();
    }

    return cpp_map;
  }

  struct PimpModel::Impl
  {
    void initialize_pimpl();

    pybind11::module_ current_module;
    pybind11::object init_f;
    pybind11::object update_f;
    pybind11::object contribution_f;
    pybind11::object division_f;
    pybind11::object show_f;
    pybind11::object get_properties_f;
    pybind11::object mass_f;
    OpaquePointer _data;
  };

  void PimpModel::Impl::initialize_pimpl()
  {
    current_module = pybind11::module_::import("modules.modules");
    init_f = current_module.attr("init");
    update_f = current_module.attr("update");
    show_f = current_module.attr("show");
    contribution_f = current_module.attr("contribution");
    get_properties_f = current_module.attr("get_properties");
    division_f = current_module.attr("division");
    mass_f = current_module.attr("mass");
  }

  PimpModel::PimpModel()
  {
    if (interpreter_init)
    {
      pimpl = new PimpModel::Impl(); // NOLINT
      pimpl->initialize_pimpl();
    }
  }

  PimpModel::PimpModel(const PimpModel& rhs)
  {
    if (rhs.pimpl != nullptr)
    {
      pimpl = new PimpModel::Impl(*rhs.pimpl); // Deep copy // NOLINT
      pimpl->initialize_pimpl();
    }
    else
    {
      pimpl = nullptr;
    }
  }

  PimpModel& PimpModel::operator=(const PimpModel& rhs)
  {
    if (this == &rhs)
    {
      return *this;
    }

    delete pimpl;

    if (rhs.pimpl != nullptr)
    {
      pimpl = new PimpModel::Impl(*rhs.pimpl); // NOLINT
      pimpl->initialize_pimpl();
    }
    else
    {
      pimpl = nullptr;
    }

    return *this;
  }

  PimpModel::PimpModel(PimpModel&& rhs) noexcept
  {
    if (rhs.pimpl != nullptr)
    {
      pimpl = rhs.pimpl;
      pimpl->initialize_pimpl();
      rhs.pimpl = nullptr;
    }
    else
    {
      pimpl = nullptr;
    }
  }

  PimpModel& PimpModel::operator=(PimpModel&& rhs) noexcept
  {
    if (this == &rhs)
    {
      return *this;
    }

    delete pimpl;

    if (rhs.pimpl != nullptr)
    {
      pimpl = rhs.pimpl;
      pimpl->initialize_pimpl();
      rhs.pimpl = nullptr;
    }
    else
    {
      pimpl = nullptr;
    }

    return *this;
  }

  PimpModel::~PimpModel()
  {
    delete pimpl;
  }

  /**
  Biological model wrap
  **/

  void PimpModel::init(MC::ParticleDataHolder& p, MC::KPRNG _rng)
  {
    if (pimpl != nullptr)
    {
      auto obj = pimpl->init_f(PYTHON_FWD(p));
      pimpl->_data.ptr = new py::object(obj); // NOLINT
    }
  }

  void PimpModel::update(const double d_t,
                         MC::ParticleDataHolder& p,
                         const LocalConcentrationView& concentration,
                         MC::KPRNG _rng)
  {
    if (pimpl != nullptr)
    {
      auto bf = pybind11::buffer_info(const_cast<double*>(concentration.data()),
                                      sizeof(double),
                                      pybind11::format_descriptor<double>::format(),
                                      1,
                                      {concentration.size()},
                                      {sizeof(double)});

      pimpl->update_f(d_t,PYTHON_CST_FWD(pimpl->_data), PYTHON_FWD(p), pybind11::array_t<double>(bf));
    }
  }

  PimpModel PimpModel::division(MC::ParticleDataHolder& p, MC::KPRNG k) noexcept
  {
    PimpModel child_model;
    if (pimpl != nullptr)
    {
      auto obj = pimpl->division_f(PYTHON_CST_FWD(pimpl->_data),PYTHON_FWD(p));
      child_model.pimpl->_data.ptr = new py::object(obj); // NOLINT
    }

    return child_model;
  }

  void PimpModel::contribution(MC::ParticleDataHolder& p, const ContributionView& contribution) noexcept
  {
    if (pimpl != nullptr)
    {
      auto access_contribs = contribution.subview();
      auto* data = access_contribs.data();
      const auto cols = access_contribs.extent(0);
      const auto rows = access_contribs.extent(1);

      py::object capsule = py::capsule(data);
      auto result = py::array_t<double_t>(
          py::buffer_info(
              data,                                         // Pointer to the data
              sizeof(double_t),                             // Size of each element
              py::format_descriptor<double_t>::format(),    // Data type format
              2,                                            // Number of dimensions
              {rows, cols},                                 // Shape of the array
              {sizeof(double_t) * cols, sizeof(double_t)}), // Strides for each dimension
          capsule);

      pimpl->contribution_f(PYTHON_CST_FWD(p), PYTHON_CST_FWD(pimpl->_data), PYTHON_FWD(result));
    }
  }

  model_properties_detail_t PimpModel::get_properties() noexcept
  {
    if (pimpl != nullptr)
    {
      return convert_dict_to_map(pimpl->get_properties_f(PYTHON_CST_FWD(pimpl->_data)));
    }
    return {{"None", 1.}};
  }

  [[nodiscard]] double PimpModel::mass() const noexcept
  {
    if (pimpl != nullptr)
    {
      auto obj= pimpl->mass_f(PYTHON_CST_FWD(pimpl->_data));
      return obj.cast<double>();
    }
    return 1.;
  }

} // namespace PythonWrap