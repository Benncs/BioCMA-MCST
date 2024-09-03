#ifndef __MC_PARTICLESHPP__
#define __MC_PARTICLESHPP__

#include <mc/particles/data_holder.hpp>
#include <mc/particles/particle_model.hpp>
#include <mc/prng/prng.hpp>
#include <utility>

namespace MC
{

  template <ParticleModel _Model> class alignas(ExecInfo::cache_line_size) BaseParticle
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
      data.update(d_t, properties, concentration, globalrng);
    }

    KOKKOS_INLINE_FUNCTION BaseParticle<_Model> division()
    {
      properties.status = CellStatus::IDLE;
      auto p = data.division(properties);
    
      return BaseParticle(properties, std::move(p));
    }

    KOKKOS_INLINE_FUNCTION void contribution(ContributionView contrib)
    {
      data.contribution(properties, contrib);
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