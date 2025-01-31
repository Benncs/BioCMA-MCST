#ifndef __SIMULATION_ALIAS_HPP__
#define __SIMULATION_ALIAS_HPP__
#include <Kokkos_Core.hpp>
#include <common/kokkos_vector.hpp>
#include <cstddef>
#include <mc/events.hpp>
#include <traits/Kokkos_IterationPatternTrait.hpp>
namespace Simulation
{

  /*Following definition are related to kokkos specifc view, those types are
   * only used during cycleprocess kernel, associated functions (get/set)
   * handle data transfer if necessary */

  template <typename ExecSpace>
  using DiagonalView = Kokkos::
      View<double*, Kokkos::LayoutLeft, ExecSpace, Kokkos::MemoryTraits<Kokkos::RandomAccess>>;

  using DiagonalViewCompute = DiagonalView<ComputeSpace>;

  template <typename Space>
  using CumulativeProbabilityView = Kokkos::
      View<double**, Kokkos::LayoutRight, Space, Kokkos::MemoryTraits<Kokkos::RandomAccess>>;

  using CumulativeProbabilityViewCompute = CumulativeProbabilityView<ComputeSpace>;

  template <typename Space>
  using NeighborsView = Kokkos::View<const std::size_t**,
                                     Kokkos::LayoutStride,
                                     Space,
                                     Kokkos::MemoryTraits<Kokkos::RandomAccess>>;

  using NeighborsViewCompute = NeighborsView<ComputeSpace>;

  using LeavingFlowIndexType = Kokkos::View<std::size_t*, Kokkos::SharedHostPinnedSpace>;
  using LeavingFlowType = Kokkos::View<double*, Kokkos::SharedHostPinnedSpace>;

} // namespace Simulation

#endif