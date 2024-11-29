#include <api/api.hpp>
#include <api/api_raw.h> 
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <memory>
#include <pybind11/cast.h>
#include <pybind11/detail/common.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/pytypes.h>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

namespace py = pybind11;

std::string wrap_repr(const wrap_c_param_t& m)
{
  char* repr = nullptr;
  repr_user_param(&m, &repr); // This function perform dynamic allocation

  std::string result(repr); // Create a std::string from the char*

  free(repr);    // NOLINT free the allocated memory
  return result; // Return the string representation
}

PYBIND11_MODULE(handle_module, m) // NOLINT (Pybind11 MACRO)
{
  // Wrapping the Handle structure
  py::class_<Handle, std::shared_ptr<Handle>>(m, "Handle");

  m.def(
      "init_handle",
      [](uint32_t n_rank, uint32_t i_rank, uint64_t id, uint32_t thread_per_proces)
      {
        auto opt = Handle::init(n_rank, i_rank, id, thread_per_proces);
        if (opt.has_value())
        {
          auto* ptr = opt.value().release();
          return std::shared_ptr<Handle>(ptr);
        }
        throw std::runtime_error("Simulation handle initialisation failed");
      },
      py::arg("n_rank"),
      py::arg("i_rank"),
      py::arg("id"),
      py::arg("thread_per_proces") = 1);

  m.def("exec", &exec);
  // m.def("apply", &apply);

  m.def("apply",
        [](std::shared_ptr<Handle>& handle, bool to_load) -> std::tuple<bool, std::string>
        {
          auto rc = handle->apply(to_load);
          bool f = static_cast<bool>(rc);
          return {f, rc.get()};
        });

  m.def("register_result_path", &register_result_path);
  m.def("register_cma_path", &register_cma_path);
  m.def("register_serde", &register_serde);
  m.def("register_parameters", &register_parameters);
  m.def("register_model_name", &register_model_name);
  m.def("make_params",
        &make_params,
        py::arg("biomass_initial_concentration"),
        py::arg("final_time"),
        py::arg("delta_time"),
        py::arg("number_particle"),
        py::arg("number_exported_result"));

  py::class_<wrap_c_param_t>(m, "UserSimulationParam")
      .def(py::init<>())
      .def_readwrite("final_time", &wrap_c_param_t::final_time)
      .def_readwrite("delta_time", &wrap_c_param_t::delta_time)
      .def_readwrite("force_override", &wrap_c_param_t::force_override)
      .def_readwrite("n_thread", &wrap_c_param_t::n_thread)
      .def_readwrite("number_exported_result", &wrap_c_param_t::number_exported_result)
      .def_readwrite("recursive", &wrap_c_param_t::recursive)
      .def_readwrite("biomass_initial_concentration",
                     &wrap_c_param_t::biomass_initial_concentration)
      .def_readwrite("number_particle", &wrap_c_param_t::number_particle)
      .def("__repr__", &wrap_repr)
      // TODO Write unittest
      .def(py::pickle(
          [](const wrap_c_param_t& p) { // __getstate__
            /* Return a tuple that fully encodes the state of the object */
            return py::make_tuple(p.final_time,
                                  p.delta_time,
                                  p.force_override,
                                  p.n_thread,
                                  p.number_exported_result,
                                  p.recursive,
                                  p.biomass_initial_concentration,
                                  p.number_particle);
          },
          [](const py::tuple& t) { // __setstate__
            constexpr std::size_t n_attributes = 8;
            if (t.size() != n_attributes)
            {
              throw std::runtime_error("Invalid state!");
            }

            /* Create a new C++ instance */
            wrap_c_param_t p{};

            // NOLINTBEGIN
            // Be cafefure using array indexing
            p.final_time = t[0].cast<double>();
            p.delta_time = t[1].cast<double>();
            p.force_override = static_cast<int>(t[2].cast<bool>());
            p.n_thread = t[3].cast<int>();
            p.number_exported_result = t[4].cast<int>();
            p.recursive = static_cast<int>(t[5].cast<bool>());
            p.biomass_initial_concentration = t[6].cast<double>();
            p.number_particle = t[7].cast<int>();
            // NOLINTEND
            return p;
          }));

  // Feed

  m.def("set_feed_constant",
        [](std::shared_ptr<Handle>& handle,
           double _f,
           std::vector<double> _target,
           std::vector<std::size_t> _position,
           std::vector<std::size_t> _species)
        { handle->set_feed_constant(_f, _target, _position, _species); });
}
