#include <simulation/pc_hydro.hpp>

namespace Simulation
{
  PreCalculatedHydroState::PreCalculatedHydroState(const FlowMatrixType &_tm)
      : transition_matrix(_tm)
  {
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
