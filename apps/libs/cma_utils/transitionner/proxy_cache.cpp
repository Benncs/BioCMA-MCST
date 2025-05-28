#include <cma_utils/cache_hydro_state.hpp>
#include <transitionner/proxy_cache.hpp>
#include <transitionner/transport.hpp>
namespace CmaUtils
{
  void ProxyPreCalculatedHydroState::set_diag_transition(std::vector<double>&& diag)
  {

    DiagonalView<HostSpace> view_host(diag.data(), diag.size());
    Kokkos::resize(state.diagonal_compute, diag.size()); // FIXME
    Kokkos::deep_copy(state.diagonal_compute, view_host);
  }

  void ProxyPreCalculatedHydroState::set_transition_matrix(
      const CmaRead::FlowMap::FlowMap_const_view_t& flows_view)
  {

    state.transition_matrix = CmaUtils::get_transition_matrix(flows_view);
    set_diag_transition(CmaUtils::get_diag_transition(state.transition_matrix));
  }

  void ProxyPreCalculatedHydroState::set_transition_matrix(FlowMatrixType&& matrix)
  {
    state.transition_matrix = matrix;
  }

  void ProxyPreCalculatedHydroState::set_cumulative_probability(
      const CmaRead::Neighbors::Neighbors_const_view_t& neighbors)
  {
    state.cumulative_probability = get_cumulative_probabilities(neighbors, state.transition_matrix);
  }

} // namespace CmaUtils