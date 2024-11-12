#include <api/api.hpp>
#include <api/api_raw.h> // Include the header with your C functions.
#include <memory>
#include <pybind11/pybind11.h>
#include <stdexcept>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>

namespace py = pybind11;

PYBIND11_MODULE(handle_module, m)
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
  m.def("apply", &apply);

  m.def("apply", &apply);
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
      .def_readwrite("number_particle", &wrap_c_param_t::number_particle);

  // Feed

  m.def("set_feed_constant",
        [](std::shared_ptr<Handle>& handle,
           double _f, std::vector<double> _target, std::vector<std::size_t> _position, std::vector<std::size_t> _species) {
            handle->set_feed_constant(_f, _target, _position, _species);
        });
}
