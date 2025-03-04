#include <Kokkos_Core.hpp>
#include <cma_utils/cache_hydro_state.hpp>
#include <common/kokkos_vector.hpp>


namespace CmaUtils
{

  PreCalculatedHydroState::PreCalculatedHydroState(const FlowMatrixType& _tm)
      : transition_matrix(_tm), diagonal_compute("diagonal_compute", 0)
  {
    Kokkos::deep_copy(diagonal_compute, 0.);
  }

  PreCalculatedHydroState::PreCalculatedHydroState() : diagonal_compute("diagonal_compute", 0)
  {
    Kokkos::deep_copy(diagonal_compute, 0.);
  }

  [[nodiscard]] const FlowMatrixType& PreCalculatedHydroState::get_transition() const
  {
    return transition_matrix;
  }

}; // namespace CmaUtils
