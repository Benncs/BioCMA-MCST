#ifndef __MC_DATA_HOLDER_PARTICLE_HPP__
#define __MC_DATA_HOLDER_PARTICLE_HPP__

#include <Kokkos_Core.hpp>
#include <cmt_common/macro_constructor_assignment.hpp>
#include <common/execinfo.hpp>

namespace MC
{
  enum class CellStatus : char
  {
    IDLE,
    DEAD,
    CYTOKINESIS,
    OUT
  };

  class alignas(ExecInfo::cache_line_size) ParticleDataHolder
  {
  public:
    KOKKOS_INLINE_FUNCTION explicit ParticleDataHolder(double _weight)
        : weight(_weight)
    {
    }
 
    KOKKOS_INLINE_FUNCTION void reset()
    {
      current_container = default_container;
      current_domain = default_domain;
      random_seed = 0;
      id = 0;
      status = default_status;
      weight = default_weight;
    }

    size_t current_container = default_container;
    size_t current_domain = default_domain;
    size_t random_seed = 0;
    uint32_t id = 0;
    CellStatus status = default_status;
    double weight = default_weight;

  private:
    static constexpr size_t default_container = 0;
    static constexpr size_t default_domain = 0;
    static constexpr CellStatus default_status = CellStatus::IDLE;
    static constexpr double default_weight = 0.;
  };
} // namespace MC

#endif