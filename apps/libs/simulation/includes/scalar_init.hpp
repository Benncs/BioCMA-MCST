#ifndef __SIMULATION_SCALAR_INIT_HPP__
#define __SIMULATION_SCALAR_INIT_HPP__

#include <memory>
#include <simulation/scalar_initializer.hpp>

namespace Simulation
{
  class ScalarSimulation;
}

namespace Simulation::impl
{
  void post_init_concentration_functor(
      bool is_two_phase_flow,
      const ScalarInitializer& scalar_init,
      const std::shared_ptr<ScalarSimulation>& liquid_scalar,
      const std::shared_ptr<ScalarSimulation>& gas_scalar = nullptr);

  void post_init_concentration_file(
      bool is_two_phase_flow,
      const ScalarInitializer& scalar_init,
      const std::shared_ptr<ScalarSimulation>& liquid_scalar,
      const std::shared_ptr<ScalarSimulation>& gas_scalar = nullptr);
} // namespace Simulation::impl

#endif
