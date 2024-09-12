#ifndef __MC_PARTICLESHPP__
#define __MC_PARTICLESHPP__

#include "common/has_serialize.hpp"
#include <mc/particles/data_holder.hpp>
#include <mc/particles/particle_model.hpp>
#include <mc/prng/prng.hpp>
#include <utility>

namespace MC
{

  template <ParticleModel _Model>
  class alignas(ExecInfo::cache_line_size) BaseParticle
  {
  public:
    using Model = _Model;

    KOKKOS_INLINE_FUNCTION explicit BaseParticle(double _weight = 0) noexcept
        : properties(_weight){};

    KOKKOS_INLINE_FUNCTION void
    clearState(MC::CellStatus _status = CellStatus::IDLE) noexcept
    {
      properties.reset();
      properties.status = _status;
    }

    KOKKOS_INLINE_FUNCTION void init(KPRNG globalrng)
    {
      data.init(properties, globalrng);
    }

    KOKKOS_INLINE_FUNCTION void
    update(double d_t,
           const LocalConcentrationView &concentration,
           KPRNG globalrng)
    {
      properties.hydraulic_time += d_t;
      properties.interdivision_time += d_t;
      data.update(d_t, properties, concentration, globalrng);
    }

    KOKKOS_INLINE_FUNCTION BaseParticle<_Model> division()
    {
      this->properties.status = CellStatus::IDLE;
      this->properties.interdivision_time = 0;
      auto p = data.division(properties);
      auto prop_child = this->properties;
      prop_child.hydraulic_time = 0;
      return BaseParticle(std::move(prop_child), std::move(p));
    }

    KOKKOS_INLINE_FUNCTION void contribution(ContributionView contrib)
    {
      data.contribution(properties, contrib);
    }

    template <class Archive> void serde(Archive &ar)
    {
      properties.serde(ar);
      if constexpr (has_serialize<_Model, Archive>())
      {
        data.serde(ar);
      }
    }

    ParticleDataHolder properties;
    _Model data{};

  private:
    BaseParticle(ParticleDataHolder props, _Model &&_model)
        : properties(props), data(std::move(_model))
    {
    }
  };

  template <ParticleModel Model> using Particle = BaseParticle<Model>;

} // namespace MC

#endif