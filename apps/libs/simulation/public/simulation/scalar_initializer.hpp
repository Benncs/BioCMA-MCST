#ifndef __SIMULATION_SCALAR_INITIALIZER_HPP__
#define __SIMULATION_SCALAR_INITIALIZER_HPP__

#include <cma_read/light_2d_view.hpp>
#include <cstddef>
#include <span>

namespace Simulation
{

  using init_scalar_f_t = void (*)(size_t, CmaRead::L2DView<double> &);
  static constexpr init_scalar_f_t default_gas_init =
      [](size_t, CmaRead::L2DView<double> &) {};

  struct ScalarInitializer
  {
    std::size_t n_species;
    std::span<double> volumesgas;
    std::span<double> volumesliq;
    init_scalar_f_t liquid_f_init;
    init_scalar_f_t gaz_f_init = default_gas_init;
    bool gas_flow = false;
  };
} // namespace Simulation

#endif