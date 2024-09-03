#ifndef __CLI_FLOWMAP_TRANSITIONNER_HPP__
#define __CLI_FLOWMAP_TRANSITIONNER_HPP__

#include <cma_read/flow_iterator.hpp>
#include <cmt_common/macro_constructor_assignment.hpp>

#include <cma_read/reactorstate.hpp>
#include <cstddef>
#include <memory>

// Foward declaration
namespace Simulation
{
  class SimulationUnit;
  class PreCalculatedHydroState;
  struct TransitionState;
} // namespace Simulation

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

    FlowMapTransitioner(
        size_t n_flowmap,
        size_t _n_per_flowmap,
        FlowmapTransitionMethod method,
        size_t number_time_step,
        std::unique_ptr<CmaRead::FlowIterator> &&iterator = nullptr,
        bool is_two_phase_flow = false);

    // FlowMapTransitioner(size_t n_flowmap,
    //                     FlowmapTransitionMethod method,
    //                     bool is_two_phase_flow = false);

    DELETE_COPY_MOVE_AC(FlowMapTransitioner)

    ~FlowMapTransitioner();

    void update_flow(Simulation::SimulationUnit &unit);
    void advance(Simulation::SimulationUnit &unit);

    void update_flow(Simulation::SimulationUnit &unit,
                     std::span<double> flows,
                     size_t n_compartment);

    [[nodiscard]] size_t get_n_timestep() const
    {
      return this->n_timestep;
    };

    [[nodiscard]] size_t getFlowIndex() const;

    // TODO REMOVE THOSE
    [[nodiscard]] size_t size() const
    {
      return iterator->size();
    }
    CmaRead::ReactorState &get_current_unchecked_mut()
    {
      return iterator->get_unchcked_mut(getFlowIndex());
    }
    [[nodiscard]] const CmaRead::ReactorState &get_current_unchecked() const
    {

      return iterator->get_unchecked(getFlowIndex());
    };

    CmaRead::ReactorState &get_unchecked_mut(size_t index)
    {
      return iterator->get_unchcked_mut(index);
    }
    [[nodiscard]] const CmaRead::ReactorState &get_unchecked(size_t index) const
    {
      return iterator->get_unchecked(index);
    };

    const CmaRead::ReactorState *getState()
    {
      return current_state;
    }

  private:
    void discontinuous_transition();
    void linear_interpolation_transition();

    bool two_phase_flow;
    std::unique_ptr<CmaRead::FlowIterator> iterator = nullptr;
    size_t n_per_flowmap;
    size_t n_flowmap;
    size_t n_timestep;
    CmaRead::ReactorState interpolated_reactor_state;

    size_t current_flowmap_count;
    size_t repetition_count;
    void calculate_full_state(const CmaRead::ReactorState &reactor_state,
                              const Simulation::SimulationUnit &unit,
                              PreCalculatedHydroState *liq_hydro_state,
                              PreCalculatedHydroState *gas_hydro_state);

    void calculate_liquid_state(
        const CmaRead::FlowMap::FlowMap_const_view_t &mat_f_liq_view,
        const Simulation::SimulationUnit &unit,
        PreCalculatedHydroState *liq_hydro_state);

    void (FlowMapTransitioner::*f_update)(Simulation::SimulationUnit &unit);

    void update_flow_interpolation(Simulation::SimulationUnit &unit);

    void update_flow_discontinous(Simulation::SimulationUnit &unit);

    std::vector<PreCalculatedHydroState> liquid_pc;
    std::vector<PreCalculatedHydroState> gas_pc;

    PreCalculatedHydroState *current_liq_hydro_state = nullptr;
    PreCalculatedHydroState *current_gas_hydro_state = nullptr;
    TransitionState *interpolated_state;

    const CmaRead::ReactorState *current_state = nullptr;

    size_t current_index;
  };

  [[nodiscard]] inline size_t FlowMapTransitioner::getFlowIndex() const
  {
    return this->repetition_count % this->n_flowmap;
  }

} // namespace Simulation

#endif //__CLI_FLOWMAP_TRANSITIONNER_HPP__
