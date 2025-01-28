#include "Kokkos_Core.hpp"
#include "common/kokkos_vector.hpp"
#include <simulation/pc_hydro.hpp>

namespace Simulation
{
  PreCalculatedHydroState::PreCalculatedHydroState(const FlowMatrixType &_tm)
      : transition_matrix(_tm)
  {

  }


  void PreCalculatedHydroState::set_diag_transition(std::vector<double>&& diag)
  {
    DiagonalView<HostSpace> view_host(diag.data(),
                             diag.size());
    view_compute = Kokkos::create_mirror_view_and_copy(ComputeSpace(),view_host);                       
  }

  // PreCalculatedHydroState &
  // PreCalculatedHydroState::operator=(PreCalculatedHydroState &&rhs) noexcept
  // {
  //   if (this != &rhs)
  //   {
  //     transition_matrix = std::move(rhs.transition_matrix);
  //   }
  //   return *this;
  // }

  // PreCalculatedHydroState::PreCalculatedHydroState(
  //     PreCalculatedHydroState &&other) noexcept
  //     : transition_matrix(std::move(other.transition_matrix))
  // {
  // }
}; // namespace Simulation
