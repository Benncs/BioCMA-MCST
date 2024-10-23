#ifndef __MC_REACTORDOMAIN_HPP__
#define __MC_REACTORDOMAIN_HPP__

#include <Kokkos_Core.hpp>
#include <Kokkos_Macros.hpp>
#include <Kokkos_StdAlgorithms.hpp>
#include <cassert>
#include <cma_read/neighbors.hpp>
#include <cmt_common/macro_constructor_assignment.hpp>
#include <common/common.hpp>
#include <common/common_types.hpp>
#include <cstddef>
#include <cstdint>
#include <mc/container_state.hpp>
#include <span>

WARN_EXPERIMENTAL

namespace MC
{

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
    SET_NON_COPYABLE(ReactorDomain)

    /**
     * @brief Default constructor
     *
     */
    ReactorDomain() = default;

    /**
     * @brief Move constructor
     **/
    ReactorDomain(ReactorDomain &&other) noexcept;

    /**
     * @brief Main constructor
     */
    ReactorDomain(std::span<double> volumes,
                  const CmaRead::Neighbors::Neighbors_const_view_t &_neighbors);
    /**
     * @brief Default destructor
     *
     */
    ~ReactorDomain() = default;

    /**
     * @brief Move assignment operator
     *
     */
    ReactorDomain &operator=(ReactorDomain &&other) noexcept;

    /**
    @brief Set volume of liquid and gas of each compartment
    */
    void setVolumes(std::span<double const> volumes_gas,
                    std::span<double const> volumes_liq);

    /**
     * @brief Update neigbors of compartments
     */
    void
    setLiquidNeighbors(const CmaRead::Neighbors::Neighbors_const_view_t &data);

    // GETTERS
    /**
     * @brief Accesses the compartment at the specified index.
     *
     * Provides non-const access to the compartment at the given index. This
     * allows modification of the compartment's state.
     */
    auto &operator[](size_t i_c);
    /**
     * @brief Accesses the compartment at the specified index.
     *
     * Provides const access to the compartment at the given index. This
     * allows modification of the compartment's state.
     */
    auto &operator[](size_t i_c) const;

    /**
     * @brief Returns an iterator to the beginning of the compartments (const
     * version).
     *
     * @return An iterator to the beginning of the compartments.
     */
    [[nodiscard]] auto begin() const;

    /**
     * @brief Returns an iterator to the end of the compartments (const
     * version).
     *
     * @return An iterator to the end of the compartments.
     */
    [[nodiscard]] auto end() const;

    /**
     * @brief Iterator to the beginning of the compartments (mutable
     * version).
     *
     * @return An iterator to the beginning of the compartments.
     */
    auto begin();

    /**
     * @brief Iterator to the end of the compartments (mutable
     * version).
     *
     * @return An iterator to the end of the compartments.
     */
    auto end();

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
    [[nodiscard]] const CmaRead::Neighbors::Neighbors_const_view_t &
    getNeighbors() const;

    /**
     * @brief Returns the number of particle per compartment
     */
    [[nodiscard]] std::vector<uint64_t> getRepartition() const;

    /**
    @brief Return a unique domain from data obtained with MPI gather
    */
    static ReactorDomain
    reduce(std::span<const size_t> data, size_t original_size, size_t n_rank);

    /**
     * @brief Get reference to raw data containers
     */
    auto &data()
    {
      return shared_containers;
    }
    // std::span<ContainerState> data();

  private:
    double _total_volume = 0.; ///< Domain total volume
    size_t id = 0;             ///< Domain ID
    size_t size = 0;           ///< Number of compartment

    /*As the containers can represent a relatively big quantity of data (number
     * of compartment up to 1000) that is access on RW either during kernel and
     * outside, containers are located into sharedSpace to limit explicit deep
     * copy.
     */
    Kokkos::View<ContainerState *, Kokkos::SharedSpace>
        shared_containers; // TODO: check with GPU if sharedspace is enough, if
                           // not use PinnedToHost, if not use explicit copy

    CmaRead::Neighbors::Neighbors_const_view_t
        neighbors; ///< Containers neighbors
  };

  inline const CmaRead::Neighbors::Neighbors_const_view_t &
  ReactorDomain::getNeighbors() const
  {

    return neighbors;
  }

  inline void ReactorDomain::setLiquidNeighbors(
      const CmaRead::Neighbors::Neighbors_const_view_t &data)
  {

    neighbors = data.to_const();
  }

  inline size_t ReactorDomain::getNumberCompartments() const noexcept
  {
    return size;
  }

  inline double ReactorDomain::getTotalVolume() const noexcept
  {
    return this->_total_volume;
  }

  inline auto &ReactorDomain::operator[](size_t i_c)
  {
    // Kokkos view is not bound checked when use release
    KOKKOS_ASSERT(i_c < size);
    return this->shared_containers(i_c);
  }
  inline auto &ReactorDomain::operator[](size_t i_c) const
  {
    KOKKOS_ASSERT(i_c < size);
    return this->shared_containers(i_c);
  }

  inline auto ReactorDomain::begin() const
  {
    return Kokkos::Experimental::begin(shared_containers);
  }
  inline auto ReactorDomain::end() const
  {
    return Kokkos::Experimental::end(shared_containers);
  }

  inline auto ReactorDomain::begin()
  {
    return Kokkos::Experimental::begin(shared_containers);
  }
  inline auto ReactorDomain::end()
  {
    return Kokkos::Experimental::end(shared_containers);
  }

} // namespace MC

#endif //__MC_REACTORDOMAIN_HPP__