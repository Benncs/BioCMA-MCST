#ifndef __MC_CONTAINER_STATE_HPP__
#define __MC_CONTAINER_STATE_HPP__

#include <common/execinfo.hpp>
#include <cstdint>
#include <Kokkos_Core.hpp>
#include <mc/particles/particle_model.hpp>

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
  struct alignas(ExecInfo::cache_line_size) ContainerState
  {
    
    double volume_liq{}; ///< Volume of liquid in the container.
    double volume_gas{}; ///< Volume of gas in the container.
    uint64_t n_cells{}; ///< Number of cells in the container. @warning Be careful
                    ///< when decrementing this value.
    uint64_t id{};    ///< Unique identifier for the container.
    // std::span<double const> concentrations; ///< Span of concentrations
                                            ///< associated with the container.
    LocalConcentrationView concentrations;
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
    template <class Archive> void serialize(Archive &ar)
    {
      ar(volume_liq, volume_gas, n_cells, id);
    }
  };
} // namespace MC

#endif //__MC_CONTAINER_STATE_HPP__