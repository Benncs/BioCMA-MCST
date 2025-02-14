#include <cma_utils/transitionner.hpp>
#include <common/common.hpp>

namespace
{
  void compute_MatFlow(const CmaRead::FlowMap::FlowMap_const_view_t& flows_view,
                       CmaUtils::PreCalculatedHydroState& matflow)
  {
    PROFILE_SECTION("host:compute_MatFlow")
    matflow.set_transition_matrix(flows_view); // FIX THIS
  }

  std::vector<double> compute_inverse_diagonal(std::span<const double> volumes)
  {
    PROFILE_SECTION("compute_inverse_diagonal")

    std::vector<double> inverse_diagonal(volumes.size());

    std::transform(volumes.begin(),
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

  void calculate_liquid_state(const CmaRead::FlowMap::FlowMap_const_view_t& mat_f_liq_view,
                              const CmaRead::Neighbors::Neighbors_const_view_t& neighbors,
                              CmaUtils::PreCalculatedHydroState& liq_hydro_state)
  {
    PROFILE_SECTION("host:calculate_liquid_state")
    compute_MatFlow(mat_f_liq_view, liq_hydro_state);
    liq_hydro_state.set_cumulative_probability(neighbors);
  }
} // namespace

namespace CmaUtils
{

  FlowMapTransitionner::FlowMapTransitionner(size_t _n_flowmap,
                                             size_t _n_per_flowmap,
                                             size_t number_time_step,
                                             std::unique_ptr<CmaRead::FlowIterator>&& _iterator,
                                             bool is_two_phase_flow)
      : two_phase_flow((is_two_phase_flow)), n_per_flowmap(_n_per_flowmap), n_flowmap(_n_flowmap),
        n_timestep(number_time_step), repetition_count(0), current_flowmap_count(0),
        iterator(std::move(_iterator))
  {
    this->liquid_pc.resize(n_flowmap);
    this->gas_pc.resize(n_flowmap);
  }

  void FlowMapTransitionner::update_flow_worker(
      std::span<double> flows,
      std::span<double> volumeLiq,
      std::span<double> volumeGas,
      const CmaRead::Neighbors::Neighbors_const_view_t& neighbors)
  {

    auto& cls = current_liq_hydro_state();

    if (need_liquid_state())
    {
      const auto mat_f_liq_view = CmaRead::FlowMap::FlowMap_const_view_t(flows, volumeLiq.size());

      calculate_liquid_state(mat_f_liq_view, neighbors, cls);
      cls.inverse_volume = compute_inverse_diagonal(volumeLiq);
      cls.volume = std::vector<double>(volumeGas.begin(), volumeGas.end());
    }
  }

  void FlowMapTransitionner::calculate_full_state(
      const CmaRead::ReactorState& reactor_state,
      CmaUtils::PreCalculatedHydroState& liq_hydro_state,
      CmaUtils::PreCalculatedHydroState& gas_hydro_state) const
  {
    calculate_liquid_state(reactor_state.liquid_flow.getViewFlows(),
                           reactor_state.liquid_flow.getViewNeighors(),
                           liq_hydro_state);

    liq_hydro_state.inverse_volume = compute_inverse_diagonal(reactor_state.liquidVolume);
    liq_hydro_state.volume = reactor_state.liquidVolume;

    if (is_two_phase_flow())
    {
      compute_MatFlow(reactor_state.gas_flow.getViewFlows(), gas_hydro_state);
      gas_hydro_state.inverse_volume = compute_inverse_diagonal(reactor_state.gasVolume);
      gas_hydro_state.volume = reactor_state.gasVolume;
    }
  }

  void FlowMapTransitionner::advance()
  {
    // unit.setLiquidFlow(current_liq_hydro_state());
    // if (iterator)
    // {
    //   if (two_phase_flow)
    //   {
    //     unit.setGasFlow(current_gas_hydro_state());
    //   }
    //   const auto& current_state = get_current_reactor_state();
    //   unit.setVolumes(current_state.gasVolume, current_state.liquidVolume);
    // }

    if (++this->current_flowmap_count == this->n_per_flowmap)
    {
      this->repetition_count++;
      this->current_flowmap_count = 0;
    }
  }

  [[nodiscard]] std::size_t FlowMapTransitionner::get_n_timestep() const noexcept
  {
    return n_timestep;
  }

} // namespace CmaUtils
