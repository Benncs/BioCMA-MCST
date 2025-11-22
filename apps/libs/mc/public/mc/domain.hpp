#ifndef __MC_REACTORDOMAIN_HPP__
#define __MC_REACTORDOMAIN_HPP__

#include <Kokkos_Core.hpp>
#include <cassert>
#include <common/common.hpp>
#include <cstddef>
#include <cstdint>
#include <mc/alias.hpp>
#include <mc/traits.hpp>
#include <span>

namespace MC
{
  template <typename Space>
  using NeighborsView =
      Kokkos::View<std::size_t**,
                   Kokkos::LayoutRight,
                   Space,
                   Kokkos::MemoryTraits<Kokkos::RandomAccess>>;

  

  /**
   * @brief Represents the spatial domain where Monte Carlo particles can exist.
   *
   * The `ReactorDomain` class models a spatial domain using a compartment-based
   * approach. Each compartment within the domain holds and manages its own
   * state through an instance of `ContainerState`. The compartments
   * collectively represent the entire spatial domain in which particles
   * interact, move, and evolve according to Monte Carlo simulations.
   *
   * The `ReactorDomain` class is responsible for maintaining the structure and
   * state of the compartments, and for providing access to the compartments for
   * simulation and modeling purposes.
   *
   * @details
   * - The domain is partitioned into multiple compartments, each represented by
   * a `ContainerState`.
   * @see ContainerState
   */
  class ReactorDomain
  {
  public:
    using n_cells_view_type = Kokkos::View<uint64_t*, ComputeSpace>;
    ReactorDomain(const ReactorDomain&) = delete;
    ReactorDomain& operator=(const ReactorDomain&) = delete;

    /**
     * @brief Default constructor
     *
     */
    ReactorDomain();

    /**
     * @brief Move constructor
     **/
    ReactorDomain(ReactorDomain&& other) noexcept = default;

    /**
     * @brief Main constructor
     */
    ReactorDomain(std::span<double> volumes, std::span<const size_t> neighbors);
    /**
     * @brief Default destructor
     *
     */
    ~ReactorDomain() = default;

    /**
     * @brief Move assignment operator
     *
     */
    ReactorDomain& operator=(ReactorDomain&& other) noexcept;

    /**
    @brief Set volume of liquid and gas of each compartment
    */
    void setVolumes(std::span<double const> volumes_liq);

    /**
     * @brief Update neigbors of compartments
     */

    void setLiquidNeighbors(std::span<const size_t> flat_data);

    /**
     * @brief Return the number of compartment in the domain
     */
    [[nodiscard]] size_t getNumberCompartments() const noexcept;

    /**
     * @brief Return total volume of domain
     */
    [[nodiscard]] double getTotalVolume() const noexcept;

    /**
     * @brief Return a const reference to neighbors
     */
    [[nodiscard]] ConstNeighborsView<ComputeSpace> getNeighbors() const;

    template <class Archive> void save(Archive& ar) const
    {

      ar(id, size, _total_volume);
    }

    template <class Archive> void load(Archive& ar)
    {
      ar(id, size, _total_volume);
    }

  private:
    double _total_volume = 0.; ///< Domain total volume
    size_t id = 0;             ///< Domain ID
    size_t size = 0;           ///< Number of compartment
    NeighborsView<ComputeSpace> k_neighbor;
  };

  inline ConstNeighborsView<ComputeSpace> ReactorDomain::getNeighbors() const
  {
    return k_neighbor;
  }

  inline size_t ReactorDomain::getNumberCompartments() const noexcept
  {
    return size;
  }

  inline double ReactorDomain::getTotalVolume() const noexcept
  {
    return this->_total_volume;
  }

} // namespace MC

#endif //__MC_REACTORDOMAIN_HPP__
