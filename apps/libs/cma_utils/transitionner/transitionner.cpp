#include <cma_utils/iteration_state.hpp>
#include <common/kokkos_vector.hpp>
#include <algorithm>
#include <transitionner/transitionner.hpp>
#include <common/common.hpp>
#include <utility>
#include <transitionner/proxy_cache.hpp>
namespace
{
  void compute_MatFlow(const CmaRead::FlowMap::FlowMap_const_view_t& flows_view,
                       CmaUtils::ProxyPreCalculatedHydroState& matflow)
  {
    PROFILE_SECTION("host:compute_MatFlow")
    matflow.set_transition_matrix(flows_view); // FIX THIS
  }

  std::vector<double> compute_inverse_diagonal(std::span<const double> volumes)
  {
    PROFILE_SECTION("compute_inverse_diagonal")

    std::vector<double> inverse_diagonal(volumes.size());

    std::ranges::transform(volumes,

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
                              CmaUtils::ProxyPreCalculatedHydroState& liq_hydro_state)
  {
    PROFILE_SECTION("host:calculate_liquid_state")
    compute_MatFlow(mat_f_liq_view, liq_hydro_state);
    liq_hydro_state.set_cumulative_probability(neighbors);
  }

} // namespace

namespace CmaUtils
{
  NeighborsView<HostSpace> FlowMapTransitionner::get_neighbors_view(
      const CmaRead::Neighbors::Neighbors_const_view_t& liquid_neighbors)
  {
    return NeighborsView<HostSpace>(const_cast<size_t*>(liquid_neighbors.data().data()),
                                    liquid_neighbors.getNRow(),
                                    liquid_neighbors.getNCol());
  }
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

    if (need_liquid_state_calculation())
    {
      auto& cls = current_liq_hydro_state();
      const auto mat_f_liq_view = CmaRead::FlowMap::FlowMap_const_view_t(flows, volumeLiq.size());
      calculate_liquid_state(mat_f_liq_view, neighbors, cls);
      cls.state.inverse_volume = compute_inverse_diagonal(volumeLiq);
      cls.state.volume = std::vector<double>(volumeLiq.begin(), volumeLiq.end());
      if (is_two_phase_flow())
      {
        auto& cgs = current_liq_hydro_state();
        cgs.state.inverse_volume = compute_inverse_diagonal(volumeGas);
        cgs.state.volume = std::vector<double>(volumeGas.begin(), volumeGas.end());
        ;
      }
    }
  }

  void FlowMapTransitionner::calculate_full_state(
      const CmaRead::ReactorState& reactor_state,
      CmaUtils::ProxyPreCalculatedHydroState& liq_hydro_state,
      CmaUtils::ProxyPreCalculatedHydroState& gas_hydro_state) const
  {
    calculate_liquid_state(reactor_state.liquid_flow.getViewFlows(),
                           reactor_state.liquid_flow.getViewNeighors(),
                           liq_hydro_state);
    // liq_hydro_state.neighbors =
    liq_hydro_state.state.inverse_volume = compute_inverse_diagonal(reactor_state.liquidVolume);
    liq_hydro_state.state.volume = reactor_state.liquidVolume;

    if (is_two_phase_flow())
    {
      compute_MatFlow(reactor_state.gas_flow.getViewFlows(), gas_hydro_state);
      gas_hydro_state.state.inverse_volume = compute_inverse_diagonal(reactor_state.gasVolume);
      gas_hydro_state.state.volume = reactor_state.gasVolume;
    }
  }

  IterationState FlowMapTransitionner::advance()
  {
    update_flow();
    auto liquid_neighbors = get_current_reactor_state().liquid_flow.getViewNeighors();
    auto host_view = NeighborsView<HostSpace>(const_cast<size_t*>(liquid_neighbors.data().data()),
                                              liquid_neighbors.getNRow(),
                                              liquid_neighbors.getNCol());

    return common_advance(host_view,{{"energy_dissipation",get_current_reactor_state().energy_dissipation}});
  }

  IterationState
  FlowMapTransitionner::advance_worker(std::span<double> flows,
                                       std::span<double> volumeLiq,
                                       std::span<double> volumeGas,
                                       const CmaRead::Neighbors::Neighbors_const_view_t& neighbors)
  {
    auto host_view = NeighborsView<HostSpace>(
        const_cast<size_t*>(neighbors.data().data()), neighbors.getNRow(), neighbors.getNCol());
    update_flow_worker(flows, volumeLiq, volumeGas, neighbors);
    
    return common_advance(host_view,{});
  }

  IterationState FlowMapTransitionner::common_advance(NeighborsView<HostSpace> host_view,std::unordered_map<std::string, std::span<const double>>&& info)
  {
    auto* const liq = &current_liq_hydro_state().state;
    auto* const gas = &current_gas_hydro_state().state;
    update_counters();
    auto& eps = get_current_reactor_state().energy_dissipation;

    

    return {.liq = liq, .gas = gas, .neighbors = std::move(host_view),.infos=info};
  }

  void FlowMapTransitionner::update_counters()
  {
    //The current index is increment
    if (++this->current_flowmap_count == this->n_per_flowmap)
    {
      //if it exceeds the number of iteration per flowmap
      //We reset it to 0 and increment repetition_count, this allow to switch to the next element in buffer 
      //During the first loop we have repetition_count<n_flowmap but after looping we have repetition_count>n_flowmap 
      this->repetition_count++;
      this->current_flowmap_count = 0;
    }
  }

  [[nodiscard]] std::size_t FlowMapTransitionner::get_n_timestep() const noexcept
  {
    return n_timestep;
  }
} // namespace CmaUtils
