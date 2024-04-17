#ifndef __CLI_UPDATE_FLOWS_HPP__
#define __CLI_UPDATE_FLOWS_HPP__

#include <cma_read/reactorstate.hpp>
#include <cstddef>
#include <simulation/simulation.hpp>
#include <simulation/vecmatflows.hpp>

namespace Simulation
{
  void update_flow(size_t &iteration_count,
                   size_t n_loop,
                   Simulation::SimulationUnit &unit,
                   ReactorState &reactor_state,
                   Simulation::BasicCacheMatflows &liquid_flows,
                   Simulation::BasicCacheMatflows &gas_flows);

  void update_flow(size_t &iteration_count,
                   size_t n_loop,
                   Simulation::SimulationUnit &unit,
                   std::span<double> flows,
                   size_t nc,
                   Simulation::BasicCacheMatflows &liq);
} // namespace Simulation

#endif //__CLI_UPDATE_FLOWS_HPP__
