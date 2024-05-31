#ifndef __CLI_UPDATE_FLOWS_HPP__
#define __CLI_UPDATE_FLOWS_HPP__

#include "birem_common/macro_constructor_assignment.hpp"
#include <cma_read/reactorstate.hpp>
#include <cstddef>
#include <simulation/basic_cache_hydro.hpp>
#include <simulation/simulation.hpp>

namespace Simulation
{

  class FlowMapTransitioner
  {
  public:
    enum FlowmapTransitionMethod
    {
      Discontinuous,
      InterpolationFO,
    };

    explicit FlowMapTransitioner(FlowmapTransitionMethod method);
    inline void perform_transition()
    {
      f_transition();
    }

  private:
    void discontinuous_transition();
    void first_order_interpolation_transition();

    std::function<void(void)> f_transition;
  };

  void update_flow(size_t &iteration_count,
                   size_t n_per_flowmap,
                   size_t n_loop,
                   Simulation::SimulationUnit &unit,
                   const ReactorState &reactor_state,
                   Simulation::BasicCacheHydro &liquid_flows,
                   Simulation::BasicCacheHydro &gas_flows,
                   bool tpf);

  void update_flow(size_t &iteration_count,
                   size_t n_loop,
                   Simulation::SimulationUnit &unit,
                   std::span<double> flows,
                   size_t nc,
                   Simulation::BasicCacheHydro &liq);
} // namespace Simulation

#endif //__CLI_UPDATE_FLOWS_HPP__
