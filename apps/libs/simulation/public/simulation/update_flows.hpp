#ifndef __CLI_UPDATE_FLOWS_HPP__
#define __CLI_UPDATE_FLOWS_HPP__

#include "birem_common/macro_constructor_assignment.hpp"
#include "cma_read/flow_iterator.hpp"
#include "simulation/pc_hydro.hpp"
#include <cma_read/reactorstate.hpp>
#include <cstddef>
#include <iomanip>
#include <memory>
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

    FlowMapTransitioner(size_t n_flowmap,
                        size_t _n_per_flowmap,
                        FlowmapTransitionMethod method,
                        size_t number_time_step,
                        std::unique_ptr<FlowIterator>&& iterator = nullptr,
                        bool is_two_phase_flow = false);

    // FlowMapTransitioner(size_t n_flowmap,
    //                     FlowmapTransitionMethod method,
    //                     bool is_two_phase_flow = false);

    DELETE_COPY_MOVE_AC(FlowMapTransitioner);

    ~FlowMapTransitioner() = default;

    void update_flow(Simulation::SimulationUnit &unit);

    void update_flow(Simulation::SimulationUnit &unit,
                     std::span<double> flows,
                     size_t n_compartment);

    [[nodiscard]] size_t get_n_timestep() const
    {
      return this->n_timestep;
    };

    // inline void perform_transition()
    // {
    //   f_transition();
    // }

    [[nodiscard]] size_t getFlowIndex() const;

    // TODO REMOVE THOSE
    [[nodiscard]] size_t size() const
    {
      return iterator->size();
    }
    ReactorState &get_current_unchecked_mut()
    {
      return iterator->get_unchcked_mut(getFlowIndex());
    }
    [[nodiscard]] const ReactorState &get_current_unchecked() const
    {
      
      return iterator->get_unchecked(getFlowIndex());
    };

    ReactorState &get_unchecked_mut(size_t index)
    {
      return iterator->get_unchcked_mut(index);
    }
    [[nodiscard]] const ReactorState &get_unchecked(size_t index) const
    {
      return iterator->get_unchecked(index);
    };

    const ReactorState* getState(){return current_state;}

  private:
    using f_transition_t = void (FlowMapTransitioner::*)();
    using f_advance_t = void (FlowMapTransitioner::*)();
    void discontinuous_transition();
    void discontinuous_advance();
    void linear_interpolation_transition();

    void linear_interpolation_advance();

    PreCalculatedHydroState linear_interpolation_pc_state(PreCalculatedHydroState& current,PreCalculatedHydroState& next);

  

    void advance(Simulation::SimulationUnit &unit);

    bool tpf;


    std::vector<TransitionState> transtion_state;
    std::unique_ptr<FlowIterator> iterator = nullptr;
    size_t n_per_flowmap;
    size_t n_flowmap;
    size_t n_timestep;

    size_t current_flowmap_count;
    size_t repetition_count;
    f_transition_t f_transition;
    void calculate_full_state(const ReactorState& reactor_state,const Simulation::SimulationUnit &unit);

    void calculate_liquid_state(const FlowMap::FlowMap_const_view_t & mat_f_liq_view,const Simulation::SimulationUnit &unit);
    f_advance_t f_advance;

    PreCalculatedHydroState *current_liq_matflow = nullptr;
    PreCalculatedHydroState *current_gas_matflow = nullptr;
    TransitionState *current_liquid_state = nullptr;
    const ReactorState* current_state =nullptr;

    
    size_t current_index;
  };

  [[nodiscard]] inline size_t FlowMapTransitioner::getFlowIndex() const
  {
    return this->repetition_count % this->n_flowmap;
  }

  

} // namespace Simulation

#endif //__CLI_UPDATE_FLOWS_HPP__
