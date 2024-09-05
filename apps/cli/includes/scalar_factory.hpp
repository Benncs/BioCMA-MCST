#ifndef __SCALAR_FACTORY__
#define __SCALAR_FACTORY__

#include <optional>
#include <simulation/scalar_initializer.hpp>
#include <span>
#include <stdexcept>
#include <variant>



namespace ScalarFactory
{
  struct Uniform
  {
    std::vector<double> liquid_concentrations;
    std::optional<std::vector<double>> gas_concentration = std::nullopt;
  };
  struct Local
  {
    std::vector<double> liquid_concentrations;
    std::vector<size_t> liquid_indices;
    std::optional<std::vector<double>> gas_concentration= std::nullopt;
    std::optional<std::vector<size_t>> gas_indices= std::nullopt;
  };
  struct File
  {
    size_t n_compartment;
    std::string_view path;
  };
  struct CustomScript
  {
    std::string_view path;
  };

  using ScalarVariant = std::variant<Uniform, Local, File, CustomScript>;

  struct Visitor
  {
    Simulation::ScalarInitializer operator()(Uniform args);
    Simulation::ScalarInitializer operator()(Local args);
    Simulation::ScalarInitializer operator()(File filepath);
    Simulation::ScalarInitializer operator()(CustomScript path);
  };

  Simulation::ScalarInitializer scalar_factory(bool f_init_gas_flow,
                                             std::span<double> gas_volume,
                                             std::span<double> liquid_volume,
                                             ScalarVariant arg_liq);

} // namespace ScalarFactory



#endif