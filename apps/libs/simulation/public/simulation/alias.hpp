#ifndef __SIMULATION_ALIAS_HPP__
#define __SIMULATION_ALIAS_HPP__
#include <Kokkos_Core.hpp>
#include <common/kokkos_vector.hpp>
#include <mc/container_state.hpp>
#include <mc/events.hpp>
#include <mc/prng/prng.hpp>
namespace Simulation
{

  /*Following definition are related to kokkos specifc view, those types are
   * only used during cycleprocess kernel, associated functions (get/set)
   * handle data transfer if necessary */
  using DiagonalViewCompute =
      Kokkos::View<double *,
                   Kokkos::LayoutLeft,
                   ComputeSpace,
                   Kokkos::MemoryTraits<Kokkos::RandomAccess>>;

  using CumulativeProbabilityViewCompute =
      Kokkos::View<double **,
                   Kokkos::LayoutLeft,
                   ComputeSpace,
                   Kokkos::MemoryTraits<Kokkos::RandomAccess>>;

  using ContributionViewCompute =
      Kokkos::View<double **, Kokkos::LayoutLeft, ComputeSpace>;
  using NeighborsViewCompute =
      Kokkos::View<const size_t **, Kokkos::LayoutStride, ComputeSpace>;

  using LeavingFlowIndexType = Kokkos::View<size_t *, ComputeSpace>;
  using LeavingFlowType = Kokkos::View<double *, ComputeSpace>;

} // namespace Simulation

#endif