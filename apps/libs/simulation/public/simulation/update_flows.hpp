#ifndef __CLI_UPDATE_FLOWS_HPP__
#define __CLI_UPDATE_FLOWS_HPP__

#include <cma_read/reactorstate.hpp>
#include <cstddef>
#include <simulation/simulation.hpp>
#include <simulation/basic_cache_hydro.hpp>

namespace Simulation
{
  void update_flow(size_t &iteration_count,
  size_t n_per_flowmap,
                   size_t n_loop,
                   Simulation::SimulationUnit &unit,
                   const ReactorState &reactor_state,
                   Simulation::BasicCacheHydro &liquid_flows,
                   Simulation::BasicCacheHydro &gas_flows,bool tpf);

  void update_flow(size_t &iteration_count,
                   size_t n_loop,
                   Simulation::SimulationUnit &unit,
                   std::span<double> flows,
                   size_t nc,
                   Simulation::BasicCacheHydro &liq);
} // namespace Simulation

#endif //__CLI_UPDATE_FLOWS_HPP__
