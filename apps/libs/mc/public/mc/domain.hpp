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

  /** @brief Store position and value of volumic flow at outlet */
  struct LeavingFlow
  {
    std::size_t index;
    double flow;
  };

  /** @brief Structure to store information about domain needed during MC cycle
    data is likely to change between each iteration
  */
  template <typename ExecSpace, bool is_const = true> struct DomainState
  {
    MC::NeighborsView<ExecSpace, is_const> neighbors;
    DiagonalView<ExecSpace, is_const> diag_transition;
    CumulativeProbabilityView<ExecSpace, is_const> cumulative_probability;
    LeavingFlowView<is_const> leaving_flow;
    VolumeView<ExecSpace, is_const> liquid_volume;
  };

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
     */
    ReactorDomain();

    /**
     * @brief Move constructor
     **/
    ReactorDomain(ReactorDomain&& other) noexcept = default;

    /**
     * @brief Main constructor
     */
    explicit ReactorDomain(std::span<double> volumes);

    ReactorDomain(double total_volume, std::size_t size);

    /**
     * @brief Default destructor
     */
    ~ReactorDomain() = default;

    /**
     * @brief Move assignment operator
     */
    ReactorDomain& operator=(ReactorDomain&& other) noexcept;

    void update(std::span<const double> newliquid_volume,
                std::span<const std::size_t> neighors_flat,
                std::span<const double> out_flows,
                std::span<const double> proba_flat);

    void set_leaving_flow(std::size_t i, std::size_t i_flow, double flow) const;

    void init_inner(std::size_t n_flows);

    [[nodiscard]] DomainState<ComputeSpace, true> get_const_inner();

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
    DomainState<ComputeSpace, false> inner;

    /**
    @brief Set volume of liquid and gas of each compartment
    */
    void setVolumes(std::span<double const> volumes_liq);

    /**
     * @brief Update neigbors of compartments
     */
    void setLiquidNeighbors(std::size_t e1,
                            std::size_t e2,
                            std::span<const size_t> flat_data);
  };

  // inline NeighborsView<ComputeSpace, true> ReactorDomain::getNeighbors() const
  // {
  //   return this->inner.neighbors;
  // }

  inline size_t ReactorDomain::getNumberCompartments() const noexcept
  {
    return size;
  }

  inline double ReactorDomain::getTotalVolume() const noexcept
  {
    return this->_total_volume;
  }

  [[nodiscard]] inline DomainState<ComputeSpace, true>
  ReactorDomain::get_const_inner()
  {
    return {inner.neighbors,
            inner.diag_transition,
            inner.cumulative_probability,
            inner.leaving_flow,
            inner.liquid_volume};
  }

} // namespace MC

#endif //__MC_REACTORDOMAIN_HPP__
