#include <cma_utils/cache_hydro_state.hpp>
#include <Kokkos_Core.hpp>
#include <common/kokkos_vector.hpp>
#include <get_cumulative_proba.hpp>
#include <cma_utils/transport.hpp>

namespace CmaUtils
{
  PreCalculatedHydroState::PreCalculatedHydroState(const FlowMatrixType& _tm)
      : transition_matrix(_tm), diagonal_compute("diagonal_compute",0)
  {
    Kokkos::deep_copy(diagonal_compute, 0.);
  }

  PreCalculatedHydroState::PreCalculatedHydroState() : diagonal_compute("diagonal_compute",0)
  {
    Kokkos::deep_copy(diagonal_compute, 0.);
  }

  void PreCalculatedHydroState::set_diag_transition(std::vector<double>&& diag)
  {

  
    DiagonalView<HostSpace> view_host(diag.data(), diag.size());
    Kokkos::resize(diagonal_compute, diag.size()); // FIXME
    Kokkos::deep_copy(diagonal_compute, view_host);
  }

  [[nodiscard]] const FlowMatrixType& PreCalculatedHydroState::get_transition() const
  {
    return transition_matrix;
  }

  void PreCalculatedHydroState::set_transition_matrix(
      const CmaRead::FlowMap::FlowMap_const_view_t& flows_view)
  {

    transition_matrix = CmaUtils::get_transition_matrix(flows_view);
    set_diag_transition(CmaUtils::get_diag_transition(transition_matrix));
  }

  void PreCalculatedHydroState::set_transition_matrix(FlowMatrixType&& matrix)
  {
    transition_matrix = matrix;
  }

  void PreCalculatedHydroState::set_cumulative_probability(
      const CmaRead::Neighbors::Neighbors_const_view_t& neighbors)
  {
    cumulative_probability = get_cumulative_probabilities(neighbors, transition_matrix);

  }

}; // namespace CmaUtils
