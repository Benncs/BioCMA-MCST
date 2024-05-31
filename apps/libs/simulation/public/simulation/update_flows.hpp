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

    explicit FlowMapTransitioner(size_t n_flowmap,size_t _n_per_flowmap,
                                 FlowmapTransitionMethod method);
    DELETE_COPY_MOVE_AC(FlowMapTransitioner);

    ~FlowMapTransitioner() = default;

    inline void perform_transition()
    {
      f_transition();
    }
    std::vector<PreCalculatedHydroState> liquid_flows;
    std::vector<PreCalculatedHydroState> gas_flows;
    size_t n_per_flowmap;
    size_t n_flowmap;

    size_t current_flownap_count;
    size_t repetition_count;

  private:
    void discontinuous_transition();
    void first_order_interpolation_transition();



    std::function<void(void)> f_transition;
  };

  void update_flow(
                   Simulation::SimulationUnit &unit,
                   const ReactorState &reactor_state,
                   FlowMapTransitioner &fmt,
                   bool tpf);

  //  void update_flow(size_t &iteration_count,
  //          size_t n_per_flowmap,
  //          size_t n_loop,
  //          Simulation::SimulationUnit &unit,
  //          const ReactorState &reactor_state,
  //          Simulation::BasicCacheHydro &liquid_flows,
  //          Simulation::BasicCacheHydro &gas_flows,
  //          bool tpf);

  void update_flow(size_t &iteration_count,
                   size_t n_loop,
                   Simulation::SimulationUnit &unit,
                   std::span<double> flows,
                   size_t n_compartment,
                   Simulation::BasicCacheHydro &liq);
} // namespace Simulation

#endif //__CLI_UPDATE_FLOWS_HPP__
