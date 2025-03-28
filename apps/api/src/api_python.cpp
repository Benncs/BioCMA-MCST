#include "core/scalar_factory.hpp"
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
#include <pybind11/stl.h>
#include <stdexcept>
#include <string>
#include <tuple>

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
  py::class_<Api::SimulationInstance, std::shared_ptr<Api::SimulationInstance>>(m, "Handle");

  m.def(
      "init_handle",
      [](const std::vector<std::string>& args)
      {
        std::vector<const char*> c_args;
        c_args.reserve(args.size());
        for (const auto& arg : args)
        {
          c_args.push_back(arg.c_str());
        }

        auto opt = Api::SimulationInstance::init(static_cast<int>(c_args.size()),
                                                 const_cast<char**>(c_args.data()));
        if (opt.has_value())
        {
          auto* ptr = opt.value().release();
          return std::shared_ptr<Api::SimulationInstance>(ptr);
        }
        throw std::runtime_error("Simulation handle initialisation failed");
      },
      py::arg("argv"));

  m.def("finalize", &finalize);
  m.def("exec", &exec);
  // m.def("apply", &apply);

  m.def("apply",
        [](std::shared_ptr<Api::SimulationInstance>& handle,
           bool to_load) -> std::tuple<bool, std::string>
        {
          auto rc = handle->apply(to_load);
          bool f = static_cast<bool>(rc);
          return {f, rc.get()};
        });

  m.def("i_rank",
        [](std::shared_ptr<Api::SimulationInstance>& handle)
        { return handle->get_exec_info().current_rank; });
  m.def("n_rank",
        [](std::shared_ptr<Api::SimulationInstance>& handle)
        { return handle->get_exec_info().n_rank; });

  m.def("register_result_path", &register_result_path);
  m.def("register_cma_path", &register_cma_path);
  m.def("register_serde", &register_serde);
  m.def("register_parameters", &register_parameters);
  m.def("register_model_name", &register_model_name);

  m.def("register_initialiser_file_path", &register_initializer_path);

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

  m.def(
      "set_liquid_feed_constant",
      [](std::shared_ptr<Api::SimulationInstance>& handle,
         double _f,
         std::vector<double> _target,
         std::vector<std::size_t> _position,
         std::vector<std::size_t> _species,
         bool fed_batch)
      { handle->set_feed_constant(_f, _target, _position, _species, false, fed_batch); },
      py::arg("handle"),
      py::arg("flow"),
      py::arg("concentration value"),
      py::arg("position"),
      py::arg("species"),
      py::arg("fed_batch") = false);

  m.def(
      "set_gas_feed_constant",
      [](std::shared_ptr<Api::SimulationInstance>& handle,
         double _f,
         std::vector<double> _target,
         std::vector<std::size_t> _position,
         std::vector<std::size_t> _species)
      { handle->set_feed_constant(_f, _target, _position, _species, true); },
      py::arg("handle"),
      py::arg("flow"),
      py::arg("concentration value"),
      py::arg("position"),
      py::arg("species"));

  m.def(
      "set_initialiser_from_data",
      [](std::shared_ptr<Api::SimulationInstance>& handle,
         std::size_t n_species,
         const py::array_t<double_t>&& py_liquid,
         std::optional<py::array_t<double_t>>&& py_gas)

      {
        auto buf = py_liquid.request();
        std::span<double> data(static_cast<double*>(buf.ptr), buf.size);
        std::vector<double> liq(data.begin(), data.end());

        std::optional<std::vector<double>> gas = std::nullopt;
        if (py_gas.has_value())
        {
          auto buf = py_gas->request();
          std::span<double> data(static_cast<double*>(buf.ptr), buf.size);
          gas = std::vector<double>(data.begin(), data.end());
        }

        handle->register_scalar_initiazer(
            Core::ScalarFactory::FullCase(n_species, std::move(liq), std::move(gas)));
      },
      py::arg("handle"),
      py::arg("n_species"),
      py::arg("liquid value"),
      py::arg("gas") = std::nullopt);
}

/**
 * @example simple_simulation.py
 *
 * This example demonstrates how to perfom a simple simulation using the python API
 *
 * @code
 *   import handle_module>
 *
 *   outfolder = "./out/"
 *   simulation_name = "my_simulation_name"
 *   cma_path = "/path/to/the/cma/" #don´t forget last /
 *   def run(params):
 *       handle = handle_module.init_simulation(outfolder,simulation_name,cma_path,params)
 *       # cma with 500 compartment, simulation with 4 species liquid only
 *       liquid_concentration_0 = np.zeros((500,4))
 *       handle_module.set_initial_concentrations(handle,liquid_concentration_0)
 *       handle_module.register_model_name(handle, "model_name")
 *        # Apply the simulation settings
 *     rc = handle_module.apply(handle, False)
 *
 *     # Check if the simulation settings were applied successfully
 *     if not rc[0]:
 *         print(rc[1])
 *         return -1
 *
 *     # Execute the simulation
 *     rc = handle_module.exec(handle);
 *
 * @endcode
 */
