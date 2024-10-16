#ifndef __MC_DATA_HOLDER_PARTICLE_HPP__
#define __MC_DATA_HOLDER_PARTICLE_HPP__

#include <Kokkos_Core.hpp>
#include <cmt_common/macro_constructor_assignment.hpp>
#include <common/execinfo.hpp>

namespace MC
{

  /**
   * @brief Enumeration that possible status for a MCParticle
   */
  enum class CellStatus : char
  {
    IDLE,        ///< No particular status (default)
    DEAD,        ///< Non active (skipped) particle
    CYTOKINESIS, ///< Particle has to divide at next cycle
    OUT          ///< Non active, moved out of domain
  };

  /**
   * @brief Main properties carried by all MonteCarlo Particle
   * These properties are handled for all particle in all cycles
   */
  class alignas(ExecInfo::cache_line_size) ParticleDataHolder
  {
  public:
    /**
     * @brief Default constructor (Can be called inside Kernel)
     */
    KOKKOS_INLINE_FUNCTION explicit ParticleDataHolder(double _weight)
        : weight(_weight)
    {
    }

    /**
     * @brief Reset to default value  (Can be called inside Kernel)
     */
    KOKKOS_INLINE_FUNCTION void reset()
    {
      current_container = default_container;
      current_domain = default_domain;
      random_seed = 0;
      id = 0;
      status = default_status;
      weight = default_weight;
    }

    size_t current_container =
        default_container; ///< Current position in the domain

    // current_domain is always 0 because current simulation only handles 1
    // domain
    size_t current_domain = default_domain; ///< In which domain particles lives

    // No used
    size_t random_seed = 0;
    uint32_t id = 0;

    CellStatus status = default_status; ///< Particle state
    double weight = default_weight;     ///< Monte-Carlo weight

  private:
    // Default values
    static constexpr size_t default_container = 0;
    static constexpr size_t default_domain = 0;
    static constexpr CellStatus default_status = CellStatus::IDLE;
    static constexpr double default_weight = 0.;
  };

} // namespace MC

#endif