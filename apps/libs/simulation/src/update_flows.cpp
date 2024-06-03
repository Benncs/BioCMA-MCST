#include "simulation/simulation.hpp"
#include <cma_read/flowmap.hpp>
#include <ctime>
#include <simulation/pc_hydro.hpp>
#include <simulation/update_flows.hpp>
#include <stdexcept>
#include <transport.hpp>

namespace Simulation
{

  static void compute_MatFlow(const FlowMap::FlowMap_const_view_t &flows_view,
                              Simulation::PreCalculatedHydroState &matflow)
  {
    const auto _mat_transition_liq =
        Simulation::get_transition_matrix(flows_view);
    matflow.transition_matrix = _mat_transition_liq;
  }

  static std::vector<double>
  compute_inverse_diagonal(std::span<const double> volumes)
  {
    std::vector<double> inverse_diagonal(volumes.size());

    std::transform(
        volumes.begin(),
        volumes.end(),
        inverse_diagonal.begin(),
        [](double volume)
        {
          if (volume == 0)
          {
            throw std::invalid_argument(
                "Setvolume: Null value of volume, matrix is not invertible");
          }
          return 1.0 / volume;
        });

    return inverse_diagonal;
  }

  FlowMapTransitioner::FlowMapTransitioner(
      size_t _n_flowmap,
      size_t _n_per_flowmap,
      FlowmapTransitionMethod method,
      size_t number_time_step,
      std::unique_ptr<FlowIterator> &&_iterator,
      bool is_two_phase_flow)
      : tpf((is_two_phase_flow)), iterator(std::move(_iterator)),
        n_per_flowmap(_n_per_flowmap), n_flowmap(_n_flowmap),
        n_timestep(number_time_step), current_flowmap_count(0),
        repetition_count(0), current_index(0)
  {

    switch (method)
    {
    case (FlowmapTransitionMethod::Discontinuous):
    {
      transtion_state.resize(n_flowmap);
      f_transition = &FlowMapTransitioner::discontinuous_transition;
      f_advance = &FlowMapTransitioner::discontinuous_advance;
      break;
    }
    case (FlowMapTransitioner::InterpolationFO):
    {
      transtion_state.resize(n_per_flowmap * n_flowmap);
      // f_transition = [this] { first_order_interpolation_transition(); };
      // f_advance = [this] { linear_interpolation_advance(); };
      break;
    }
    }
  }

  void FlowMapTransitioner::linear_interpolation_advance()
  {
    discontinuous_advance();
  }

  PreCalculatedHydroState FlowMapTransitioner::linear_interpolation_pc_state(
      PreCalculatedHydroState &current, PreCalculatedHydroState &next)
  {
    return std::move(current);
  }

  void FlowMapTransitioner::linear_interpolation_transition()
  {
    throw std::runtime_error("not implemented yet");
  }

  void FlowMapTransitioner::discontinuous_transition()
  {
    current_liq_matflow = &this->transtion_state[current_index].liquid_pc;
    if (tpf)
    {
      current_gas_matflow = &this->transtion_state[current_index].gas_pc;
    }
  }

  // ok dont modify
  void FlowMapTransitioner::discontinuous_advance()
  {
    if (++this->current_flowmap_count == this->n_per_flowmap)
    {
      this->repetition_count++;
      this->current_flowmap_count = 0;
    }
  }

  // MPI Host
  void FlowMapTransitioner::update_flow(Simulation::SimulationUnit &unit)
  {
    current_index = getFlowIndex();
    current_liq_matflow = &this->transtion_state[current_index].liquid_pc;
    if (tpf)
    {
      current_gas_matflow = &this->transtion_state[current_index].gas_pc;
    }

    if (this->repetition_count < this->n_flowmap)
    {
      this->transtion_state[current_index].state =
          &get_unchecked(current_index);
      current_state = this->transtion_state[current_index].state;
      calculate_full_state(*current_state, unit);
    }
    current_state = this->transtion_state[current_index].state;

    unit.mc_unit->domain.setLiquidNeighbors(
        current_state->liquid_flow.getViewNeighors());

    advance(unit);
  }

  // MPI worker
  void FlowMapTransitioner::update_flow(Simulation::SimulationUnit &unit,
                                        std::span<double> flows,
                                        size_t n_compartment)
  {
    (this->*f_transition)();

    if (repetition_count < n_flowmap)
    {
      const auto mat_f_liq_view =
          FlowMap::FlowMap_const_view_t(flows, n_compartment);
      calculate_liquid_state(mat_f_liq_view, unit);
    }

    advance(unit);
  }

  // ok dont modify
  void FlowMapTransitioner::advance(Simulation::SimulationUnit &unit)
  {

    unit.setLiquidFlow(current_liq_matflow);
    if (iterator)
    {
      if (tpf)
      {
        unit.setGasFlow(current_gas_matflow);
      }

      unit.setVolumes(current_state->gasVolume, current_state->liquidVolume);
    }

    (this->*f_advance)();
  }

  // ok dont modify
  void FlowMapTransitioner::calculate_liquid_state(
      const FlowMap::FlowMap_const_view_t &mat_f_liq_view,
      const Simulation::SimulationUnit &unit)
  {
    compute_MatFlow(mat_f_liq_view, *current_liq_matflow);
    current_liq_matflow->cumulative_probability =
        get_CP(unit.mc_unit->domain.getNeighbors(),
               current_liq_matflow->transition_matrix);
  }
  // ok dont modify
  void FlowMapTransitioner::calculate_full_state(
      const ReactorState &reactor_state, const Simulation::SimulationUnit &unit)
  {
    calculate_liquid_state(reactor_state.liquid_flow.getViewFlows(), unit);

    current_liq_matflow->inverse_volume =
        compute_inverse_diagonal(reactor_state.liquidVolume);

    if (tpf)
    {
      compute_MatFlow(reactor_state.gas_flow.getViewFlows(),
                      *current_gas_matflow);
      current_gas_matflow->inverse_volume =
          compute_inverse_diagonal(reactor_state.gasVolume);
    }
  }

} // namespace Simulation