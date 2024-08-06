#include "cma_read/flowinfo.hpp"
#include "cma_read/reactorstate.hpp"
#include "simulation/simulation.hpp"
#include <cma_read/flowmap.hpp>
#include <cmath>
#include <ctime>
#include <simulation/pc_hydro.hpp>
#include <simulation/update_flows.hpp>
#include <stdexcept>
#include <transport.hpp>
#include <get_cumulative_proba.hpp>

double linter(double a, double b, double t)
{
  return (1 - t) * a + t * b;
}

namespace Simulation
{

  static void
  compute_MatFlow(const CmaRead::FlowMap::FlowMap_const_view_t &flows_view,
                  Simulation::PreCalculatedHydroState &matflow)
  {
    const auto _mat_transition_liq =
        Simulation::get_transition_matrix(flows_view);
    matflow.transition_matrix = _mat_transition_liq;

    matflow.diag_transition = get_diag_transition(matflow.transition_matrix);
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
      std::unique_ptr<CmaRead::FlowIterator> &&_iterator,
      bool is_two_phase_flow)
      : two_phase_flow((is_two_phase_flow)), iterator(std::move(_iterator)),
        n_per_flowmap(_n_per_flowmap), n_flowmap(_n_flowmap),
        n_timestep(number_time_step), current_flowmap_count(0),
        repetition_count(0), current_index(0)
  {
    // transtion_state.resize(n_flowmap);
    this->liquid_pc.resize(n_flowmap);
    this->gas_pc.resize(n_flowmap);
    switch (method)
    {
    case (FlowmapTransitionMethod::Discontinuous):
    {
      std::cout << "Using Discontinous transition\r\n";
      f_update = &FlowMapTransitioner::update_flow_discontinous;
      break;
    }
    case (FlowMapTransitioner::InterpolationFO):
    {
      std::cout << "Using linear interpolation transition\r\n";
      f_update = &FlowMapTransitioner::update_flow_interpolation;
      break;
    }
    }
  }

  void linter_eigen(const PreCalculatedHydroState &current,
                    const PreCalculatedHydroState &next,
                    PreCalculatedHydroState &interpolated,
                    double t)
  {
    interpolated.cumulative_probability =
        (1 - t) * current.cumulative_probability +
        t * next.cumulative_probability;

    interpolated.transition_matrix =
        (1.0 - t) * current.transition_matrix + t * next.transition_matrix;
  }

  void FlowMapTransitioner::linear_interpolation_transition()
  {
    size_t next_index = (this->repetition_count + 1) % this->n_flowmap;
    double t = static_cast<double>(this->current_flowmap_count) /
               static_cast<double>(this->n_per_flowmap);

    auto &current_l_state = this->liquid_pc[current_index];
    auto &next_l_state = this->liquid_pc[next_index];

    linter_eigen(
        current_l_state, next_l_state, interpolated_state.liquid_pc, t);

    if (two_phase_flow)
    {
      auto &current_g_state = this->gas_pc[current_index];
      auto &next_g_state = this->gas_pc[next_index];

      linter_eigen(current_g_state, next_g_state, interpolated_state.gas_pc, t);
    }

    auto n_compartment = current_state->n_compartments;
    // TODO CHECK move assigment
    interpolated_state.state.liquid_flow = current_state->liquid_flow;
    const auto *current_r = &get_unchecked(current_index);
    const auto *next_r = &get_unchecked(next_index);

    interpolated_state.state.liquidVolume.resize(n_compartment);
    interpolated_state.state.gasVolume.resize(n_compartment);
    interpolated_state.state.energy_dissipation.resize(n_compartment);
    interpolated_state.liquid_pc.inverse_volume.resize(
        current_l_state.inverse_volume.size());

    for (size_t is = 0; is < n_compartment; ++is)
    {

      interpolated_state.liquid_pc.inverse_volume[is] =
          linter(current_l_state.inverse_volume[is],
                 next_l_state.inverse_volume[is],
                 t);

      if (two_phase_flow)
      {
        interpolated_state.gas_pc.inverse_volume[is] =
            linter(this->gas_pc[current_index].inverse_volume[is],
                   this->gas_pc[next_index].inverse_volume[is],
                   t);
      }

      interpolated_state.state.liquidVolume[is] =
          linter(current_r->liquidVolume[is], next_r->liquidVolume[is], t);
      interpolated_state.state.gasVolume[is] =
          linter(current_r->gasVolume[is], next_r->gasVolume[is], t);
      interpolated_state.state.energy_dissipation[is] = linter(
          current_r->energy_dissipation[is], next_r->energy_dissipation[is], t);
    }

    interpolated_state.state.n_compartments = n_compartment;
    current_liq_hydro_state = &interpolated_state.liquid_pc;
    current_gas_hydro_state = &interpolated_state.gas_pc;
    current_state = &(interpolated_state.state);
  }

  void FlowMapTransitioner::update_flow_interpolation(
      Simulation::SimulationUnit &unit)
  {
    current_index = getFlowIndex();
    if (this->repetition_count < this->n_flowmap - 1)
    {

      const auto *tmp_current_sate = &get_unchecked(current_index);
      unit.mc_unit->domain.setLiquidNeighbors(
          tmp_current_sate->liquid_flow.getViewNeighors());
      calculate_full_state(*tmp_current_sate,
                           unit,
                           &this->liquid_pc[current_index],
                           &this->gas_pc[current_index]);

      // Calculate next state
      auto next_index = current_index + 1;
      const auto *next_state = &get_unchecked(next_index);
      unit.mc_unit->domain.setLiquidNeighbors(
          next_state->liquid_flow.getViewNeighors());
      calculate_full_state(*next_state,
                           unit,
                           &this->liquid_pc[next_index],
                           &this->gas_pc[next_index]);

      current_state = tmp_current_sate;
    }

    linear_interpolation_transition();
    unit.mc_unit->domain.setLiquidNeighbors(
        current_state->liquid_flow.getViewNeighors());
  }

  void FlowMapTransitioner::update_flow_discontinous(
      Simulation::SimulationUnit &unit)
  {
    current_index = getFlowIndex();

    current_liq_hydro_state = &this->liquid_pc[current_index];
    if (two_phase_flow)
    {
      current_gas_hydro_state = &this->gas_pc[current_index];
    }

    if (this->repetition_count < this->n_flowmap)
    {
      current_state = &get_unchecked(current_index);

      calculate_full_state(*current_state,
                           unit,
                           current_liq_hydro_state,
                           current_gas_hydro_state);
    }

    current_state = &get_unchecked(current_index);

    unit.mc_unit->domain.setLiquidNeighbors(
        current_state->liquid_flow.getViewNeighors());
  }

  // MPI Host
  void FlowMapTransitioner::update_flow(Simulation::SimulationUnit &unit)
  {
    (this->*f_update)(unit);
  }

  // MPI worker
  void FlowMapTransitioner::update_flow(Simulation::SimulationUnit &unit,
                                        std::span<double> flows,
                                        size_t n_compartment)
  {
    current_index = getFlowIndex();

    current_liq_hydro_state = &this->liquid_pc[current_index];
    if (two_phase_flow)
    {
      current_gas_hydro_state = &this->gas_pc[current_index];
    }

    if (repetition_count < n_flowmap)
    {
      const auto mat_f_liq_view =
          CmaRead::FlowMap::FlowMap_const_view_t(flows, n_compartment);
      calculate_liquid_state(mat_f_liq_view, unit, current_liq_hydro_state);
    }
  }

  // ok dont modify
  void FlowMapTransitioner::advance(Simulation::SimulationUnit &unit)
  {

    unit.setLiquidFlow(current_liq_hydro_state);
    if (iterator)
    {
      if (two_phase_flow)
      {
        unit.setGasFlow(current_gas_hydro_state);
      }

      unit.setVolumes(current_state->gasVolume, current_state->liquidVolume);
    }

    if (++this->current_flowmap_count == this->n_per_flowmap)
    {
      this->repetition_count++;
      this->current_flowmap_count = 0;
    }
  }

  // ok dont modify
  void FlowMapTransitioner::calculate_liquid_state(
      const CmaRead::FlowMap::FlowMap_const_view_t &mat_f_liq_view,
      const Simulation::SimulationUnit &unit,
      PreCalculatedHydroState *liq_hydro_state)
  {
    compute_MatFlow(mat_f_liq_view, *liq_hydro_state);
    liq_hydro_state->cumulative_probability =
        get_cumulative_probabilities(unit.mc_unit->domain.getNeighbors(),
               liq_hydro_state->transition_matrix);
  }
  // ok dont modify
  void FlowMapTransitioner::calculate_full_state(
      const CmaRead::ReactorState &reactor_state,
      const Simulation::SimulationUnit &unit,
      PreCalculatedHydroState *liq_hydro_state,
      PreCalculatedHydroState *gas_hydro_state)
  {
    calculate_liquid_state(
        reactor_state.liquid_flow.getViewFlows(), unit, liq_hydro_state);

    liq_hydro_state->inverse_volume =
        compute_inverse_diagonal(reactor_state.liquidVolume);

    if (two_phase_flow)
    {
      compute_MatFlow(reactor_state.gas_flow.getViewFlows(), *gas_hydro_state);
      gas_hydro_state->inverse_volume =
          compute_inverse_diagonal(reactor_state.gasVolume);
    }
  }

} // namespace Simulation