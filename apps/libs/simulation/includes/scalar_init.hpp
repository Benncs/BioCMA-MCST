#ifndef __SIMULATION_SCALAR_INIT_HPP__
#define __SIMULATION_SCALAR_INIT_HPP__

#include <memory>
#include <simulation/scalar_initializer.hpp>

namespace Simulation::impl
{
  class ScalarSimulation;
}

namespace Simulation::impl
{
  void post_init_concentration_functor(
      bool is_two_phase_flow,
      ScalarInitializer&& scalar_init,
      const std::shared_ptr<ScalarSimulation>& scalar);

  void
  post_init_concentration_file(bool is_two_phase_flow,
                               ScalarInitializer&& scalar_init,
                               const std::shared_ptr<ScalarSimulation>& scalar,
                               bool change_layout = true);
} // namespace Simulation::impl

#endif
