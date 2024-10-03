#ifndef __SIMULATION_ALIAS_HPP__
#define __SIMULATION_ALIAS_HPP__
#include <Kokkos_Core.hpp>
#include <Kokkos_Core_fwd.hpp>
#include <common/kokkos_vector.hpp>
#include <mc/container_state.hpp>
#include <mc/events.hpp>
#include <mc/prng/prng.hpp>
#include <traits/Kokkos_IterationPatternTrait.hpp>
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
                   Kokkos::LayoutRight,
                   ComputeSpace,
                   Kokkos::MemoryTraits<Kokkos::RandomAccess>>;

  //   using ContributionViewCompute =
  //       Kokkos::View<double **,
  //                    Kokkos::LayoutLeft,
  //                    ComputeSpace,
  //                    Kokkos::MemoryTraits<Kokkos::RandomAccess>>;
  // using ContributionViewCompute =
  //   Kokkos::Experimental::ScatterView<double **,
  //                Kokkos::LayoutLeft,
  //                ComputeSpace,
  //                Kokkos::MemoryTraits<Kokkos::RandomAccess>>;

  using NeighborsViewCompute =
      Kokkos::View<const size_t **,
                   Kokkos::LayoutRight,
                   ComputeSpace,
                   Kokkos::MemoryTraits<Kokkos::RandomAccess>>;

  using LeavingFlowIndexType =
      Kokkos::View<size_t *, Kokkos::SharedHostPinnedSpace>;
  using LeavingFlowType = Kokkos::View<double *, Kokkos::SharedHostPinnedSpace>;

} // namespace Simulation

#endif