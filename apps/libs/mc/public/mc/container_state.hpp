#ifndef __MC_CONTAINER_STATE_HPP__
#define __MC_CONTAINER_STATE_HPP__

#include <Kokkos_Core.hpp>
#include <common/execinfo.hpp>
#include <cstdint>
#include <mc/particles/particle_model.hpp>

using LocalConcentrationView = Kokkos::Subview<Kokkos::View<const double**>,
                                               int,
                                               decltype(Kokkos::ALL)>; ///< Concentration inside the
                                                                       ///< current container that
                                                                       ///< particle can access //TODO REMOVE 


namespace MC
{
  /**
   * @struct ContainerState
   * @brief A data structure representing the state of a container in the
   * simulation.
   *
   * The `ContainerState` struct stores information about the liquid and gas
   * volumes, the number of cells, a unique identifier, and a span of
   * concentrations. It is aligned to the cache line size specified by
   * `ExecInfo::cache_line_size`.
   */
  struct ContainerState
  {
    uint64_t id{};       ///< Unique identifier for the container.
    LocalConcentrationView concentrations;
    // alignas(ExecInfo::cache_line_size)
    //     int64_t n_cells{}; ///< Number of cells in the container. @warning Be careful
    //                        ///< when decrementing this value.

    /**
     * @brief Serializes the `ContainerState` data members.
     *
     * This function is used to serialize the state of the container, excluding
     * the `concentrations` span, which is not serialized due to its nature as a
     * view.
     *
     * @tparam Archive The type of the archive used for serialization.
     * @param ar Reference to the archive object.
     */
    template <class Archive> void serialize(Archive& ar)
    {
      ar(id);
    }

    ContainerState() = default;
  };

  static_assert(sizeof(ContainerState) <= 2 * ExecInfo::cache_line_size, "Container State size");

} // namespace MC

#endif //__MC_CONTAINER_STATE_HPP__
