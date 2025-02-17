#ifndef __SIMULATION_SCALAR_INITIALIZER_HPP__
#define __SIMULATION_SCALAR_INITIALIZER_HPP__

#include <cstddef>
#include <cstdint>
#include <functional>
#include <optional>
#include <span>
#include <vector> 

namespace Simulation
{
  enum class ScalarInitialiserType:std::uint8_t
  {
    Uniform,     // Need vec{uniform concentration}
    Local,       // Need vec{index compartment }+ vec{concentration}
    File,        // Need filepath
    CustomScript // Need script path
  };

  // using init_scalar_f_t =
  //     std::function<void(size_t, CmaRead::L2DView<double> &)>;

  using init_scalar_f_t =
      std::function<double(std::size_t, std::size_t)>;


  struct ScalarInitializer
  {
    std::size_t n_species;
    std::span<double> volumesgas;
    std::span<double> volumesliq;
    ScalarInitialiserType type;

    std::optional<init_scalar_f_t> liquid_f_init = std::nullopt;
    std::optional<init_scalar_f_t> gas_f_init = std::nullopt;
    std::optional<std::vector<double>> liquid_buffer = std::nullopt;
    std::optional<std::vector<double>> gas_buffer = std::nullopt;
    bool gas_flow = false;
  };

} // namespace Simulation

#endif