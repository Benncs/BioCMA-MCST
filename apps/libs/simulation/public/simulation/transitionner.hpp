#ifndef __CLI_FLOWMAP_TRANSITIONNER_HPP__
#define __CLI_FLOWMAP_TRANSITIONNER_HPP__

#include <cma_read/flow_iterator.hpp>
#include <cma_read/reactorstate.hpp>
#include <cstddef>
#include <memory>

namespace CmaUtils
{
  class PreCalculatedHydroState;
  
} // namespace CmaUtils

// Foward declaration
namespace Simulation
{
  class SimulationUnit;
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

    FlowMapTransitioner(size_t n_flowmap,
                        size_t _n_per_flowmap,
                        FlowmapTransitionMethod method,
                        size_t number_time_step,
                        std::unique_ptr<CmaRead::FlowIterator>&& iterator = nullptr,
                        bool is_two_phase_flow = false);

    FlowMapTransitioner(const FlowMapTransitioner&) = delete;
    FlowMapTransitioner(FlowMapTransitioner&&) = delete;
    FlowMapTransitioner& operator=(const FlowMapTransitioner&) = delete;
    FlowMapTransitioner& operator=(FlowMapTransitioner&&) = delete;

    ~FlowMapTransitioner();

    void update_flow(Simulation::SimulationUnit& unit);
    void advance(Simulation::SimulationUnit& unit);

    void update_flow(std::span<double> flows,
                     size_t n_compartment,
                     const CmaRead::Neighbors::Neighbors_const_view_t& neighbors);

    [[nodiscard]] size_t get_n_timestep() const;

    [[nodiscard]] size_t getFlowIndex() const;

    // TODO REMOVE THOSE
    [[nodiscard]] size_t size() const;

    CmaRead::ReactorState& get_current_unchecked_mut();

    [[nodiscard]] const CmaRead::ReactorState& get_current_unchecked() const;

    CmaRead::ReactorState& get_unchecked_mut(size_t index);
    [[nodiscard]] const CmaRead::ReactorState& get_unchecked(size_t index) const;

    [[nodiscard]] const CmaRead::ReactorState* getState() const;

  private:
    void discontinuous_transition();
    void linear_interpolation_transition();

    void calculate_full_state(const CmaRead::ReactorState& reactor_state,
                              const Simulation::SimulationUnit& unit,
                              CmaUtils::PreCalculatedHydroState* liq_hydro_state,
                              CmaUtils::PreCalculatedHydroState* gas_hydro_state);

    void calculate_liquid_state(const CmaRead::FlowMap::FlowMap_const_view_t& mat_f_liq_view,
                                const CmaRead::Neighbors::Neighbors_const_view_t& neighbors,
                                CmaUtils::PreCalculatedHydroState* liq_hydro_state);

    void (FlowMapTransitioner::*f_update)(Simulation::SimulationUnit& unit);

    void update_flow_interpolation(Simulation::SimulationUnit& unit);

    void update_flow_discontinous(Simulation::SimulationUnit& unit);

    bool two_phase_flow;
    std::unique_ptr<CmaRead::FlowIterator> iterator = nullptr;
    size_t n_per_flowmap;
    size_t n_flowmap;
    size_t n_timestep;
    CmaRead::ReactorState interpolated_reactor_state;

    size_t current_flowmap_count;
    size_t repetition_count;

    std::vector<CmaUtils::PreCalculatedHydroState> liquid_pc;
    std::vector<CmaUtils::PreCalculatedHydroState> gas_pc;

    CmaUtils::PreCalculatedHydroState* current_liq_hydro_state = nullptr;
    CmaUtils::PreCalculatedHydroState* current_gas_hydro_state = nullptr;
    TransitionState* interpolated_state;

    const CmaRead::ReactorState* current_state = nullptr;
    size_t current_index;
  };

  inline const CmaRead::ReactorState* FlowMapTransitioner::getState() const
  {
    return current_state;
  }

  inline CmaRead::ReactorState& FlowMapTransitioner::get_unchecked_mut(size_t index)
  {
    return iterator->get_unchcked_mut(index);
  }
  [[nodiscard]] inline const CmaRead::ReactorState&
  FlowMapTransitioner::get_unchecked(size_t index) const
  {
    return iterator->get_unchecked(index);
  };

  [[nodiscard]] inline const CmaRead::ReactorState&
  FlowMapTransitioner::get_current_unchecked() const
  {
    return iterator->get_unchecked(getFlowIndex());
  };

  inline CmaRead::ReactorState& FlowMapTransitioner::get_current_unchecked_mut()
  {
    return iterator->get_unchcked_mut(getFlowIndex());
  }

  [[nodiscard]] inline size_t FlowMapTransitioner::get_n_timestep() const
  {
    return this->n_timestep;
  };

  [[nodiscard]] inline size_t FlowMapTransitioner::getFlowIndex() const
  {
    return this->repetition_count % this->n_flowmap;
  }

  [[nodiscard]] inline size_t FlowMapTransitioner::size() const
  {
    return iterator->size();
  }

} // namespace Simulation

#endif //__CLI_FLOWMAP_TRANSITIONNER_HPP__
