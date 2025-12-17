#ifndef __SIMULATION_ALIAS_HPP__
#define __SIMULATION_ALIAS_HPP__
#include "mc/alias.hpp"
#include <Kokkos_Core.hpp>

#include <cstddef>
namespace Simulation
{

  /*Following definition are related to kokkos specifc view, those types are
   * only used during cycleprocess kernel, associated functions (get/set)
   * handle data transfer if necessary */

  // template <typename ExecSpace>
  // using DiagonalView = MC::DiagonalView<ExecSpace>;

  // template <typename ExecSpace>
  // using CumulativeProbabilityView = MC::CumulativeProbabilityView<ExecSpace>;

  using LeavingFlowIndexType
      = Kokkos::View<std::size_t*, Kokkos::SharedHostPinnedSpace>;
  using LeavingFlowType = Kokkos::View<double*, Kokkos::SharedHostPinnedSpace>;

} // namespace Simulation

#endif