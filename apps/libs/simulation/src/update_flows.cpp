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
      std::unique_ptr<FlowIterator>&& _iterator,
      bool is_two_phase_flow)
      : tpf((is_two_phase_flow)), iterator(std::move(_iterator)),
        n_per_flowmap(_n_per_flowmap), n_flowmap(_n_flowmap),
        n_timestep(number_time_step), current_flowmap_count(0),
        repetition_count(0)
  {
    liquid_flows.resize(n_flowmap);
    gas_flows.resize(n_flowmap);
    
    switch (method)
    {
    case (FlowmapTransitionMethod::Discontinuous):
    {
      liquid_transtion_state.resize(n_flowmap);
      f_transition = [this] { discontinuous_transition(); };
      f_advance = [this] { discontinuous_advance(); };
      break;
    }
    case (FlowMapTransitioner::InterpolationFO):
    {
      liquid_transtion_state.resize(n_per_flowmap * n_flowmap);
      f_transition = [this] { first_order_interpolation_transition(); };
      break;
    }
    }
  }

  void FlowMapTransitioner::first_order_interpolation_transition()
  {
    throw std::runtime_error("not implemented yet");
  }

  void FlowMapTransitioner::discontinuous_transition()
  {
    size_t current_index_mat = this->repetition_count % this->n_flowmap;
    current_liq_matflow = &this->liquid_flows[current_index_mat];
    if (tpf)
    {
      current_gas_matflow = &this->gas_flows[current_index_mat];
    }
  }

  void FlowMapTransitioner::discontinuous_advance()
  {
    if (++this->current_flowmap_count == this->n_per_flowmap)
    {
      this->repetition_count++;
      this->current_flowmap_count = 0;
    }
  }

  void FlowMapTransitioner::update_flow(Simulation::SimulationUnit &unit,
                                        std::span<double> flows,
                                        size_t n_compartment)
  {
    f_transition();

    if (repetition_count < n_flowmap)
    {
      const auto mat_f_liq_view =
          FlowMap::FlowMap_const_view_t(flows, n_compartment);
      compute_MatFlow(mat_f_liq_view, *current_liq_matflow);
      current_liq_matflow->cumulative_probability =
          get_CP(unit.mc_unit->domain.getNeighbors(),
                 current_liq_matflow->transition_matrix);
    }

    advance(unit);
  }

  void FlowMapTransitioner::update_flow(Simulation::SimulationUnit &unit,
                                        const ReactorState &reactor_state)
  {
    f_transition();

    unit.mc_unit->domain.setLiquidNeighbors(
        reactor_state.liquid_flow.getViewNeighors());

    if (this->repetition_count < this->n_flowmap)
    {

      compute_MatFlow(reactor_state.liquid_flow.getViewFlows(),
                      *current_liq_matflow);

      current_liq_matflow->inverse_volume =
          compute_inverse_diagonal(reactor_state.liquidVolume);

      if (tpf)
      {
        compute_MatFlow(reactor_state.gas_flow.getViewFlows(),
                        *current_gas_matflow);
        current_gas_matflow->inverse_volume =
            compute_inverse_diagonal(reactor_state.gasVolume);
      }

      current_liq_matflow->cumulative_probability =
          get_CP(unit.mc_unit->domain.getNeighbors(),
                 current_liq_matflow->transition_matrix);
    }

    // auto gv = reactor_state.gasVolume;
    // auto lv = reactor_state.liquidVolume;
    //  unit.setVolumes(gv,lv);
    advance(unit);
  }

  void FlowMapTransitioner::advance(Simulation::SimulationUnit &unit)
  {

    unit.setLiquidFlow(current_liq_matflow);
    if (tpf)
    {
      unit.setGasFlow(current_gas_matflow);
    }

    f_advance();
  }

} // namespace Simulation