#ifndef __CMA_UTILS_HPP__
#define __CMA_UTILS_HPP__

#include <cma_utils/cache_hydro_state.hpp>
#include <string>
#include <unordered_map>

namespace CmaUtils
{
  template <typename Space>
  using NeighborsView = Kokkos::
      View<std::size_t**, Kokkos::LayoutRight, Space, Kokkos::MemoryTraits<Kokkos::RandomAccess>>;

  /**
   * @brief Structure to store information about the reactor state during simulation.
   *
   * This structure holds data that is not directly calculated by the program but is derived
   * from compartment properties.
   */
  struct IterationState
  {
    CmaUtils::PreCalculatedHydroState* liq; ///< Pointer to the liquid phase hydrodynamic state.
    CmaUtils::PreCalculatedHydroState* gas; ///< Pointer to the gas phase hydrodynamic state.
    NeighborsView<HostSpace> neighbors;     ///< View of neighboring compartments.
    std::unordered_map<std::string, std::span<const double>>
        infos; ///< Additional information mapped by string keys.
  };
} // namespace CmaUtils

#endif