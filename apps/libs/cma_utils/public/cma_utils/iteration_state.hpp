#ifndef __CMA_UTILS_HPP__
#define __CMA_UTILS_HPP__

#include <cma_utils/cache_hydro_state.hpp>

namespace CmaUtils
{
  template <typename Space>
  using NeighborsView = Kokkos::
      View<std::size_t**, Kokkos::LayoutRight, Space, Kokkos::MemoryTraits<Kokkos::RandomAccess>>;
  struct IterationState
  {
    CmaUtils::PreCalculatedHydroState* liq;
    CmaUtils::PreCalculatedHydroState* gas;
    NeighborsView<HostSpace> neighbors;
  };
} // namespace CmaUtils

#endif