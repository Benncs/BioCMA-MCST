#include "Kokkos_Core.hpp"
#include "common/kokkos_vector.hpp"
#include <pc_hydro.hpp>
#include <transport.hpp>

namespace Simulation
{
  PreCalculatedHydroState::PreCalculatedHydroState(const FlowMatrixType& _tm)
      : transition_matrix(_tm)
  {
  }

  void PreCalculatedHydroState::set_diag_transition(std::vector<double>&& diag)
  {
    DiagonalView<HostSpace> view_host(diag.data(), diag.size());
    diagonal_compute = Kokkos::create_mirror_view_and_copy(ComputeSpace(), view_host);
  }

  [[nodiscard]] const FlowMatrixType& PreCalculatedHydroState::get_transition() const
  {
    return transition_matrix;
  }

  void PreCalculatedHydroState::set_transition_matrix(
      const CmaRead::FlowMap::FlowMap_const_view_t& flows_view)
  {

    transition_matrix = Simulation::get_transition_matrix(flows_view);
    set_diag_transition(Simulation::get_diag_transition(transition_matrix));
  }

}; // namespace Simulation
