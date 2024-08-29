#ifndef __MC_PARTICLESHPP__
#define __MC_PARTICLESHPP__

#include "common/kokkos_vector.hpp"
#include "mc/prng/prng.hpp"
#include <Kokkos_Printf.hpp>
#include <mc/particles/data_holder.hpp>
#include <mc/particles/particle_model.hpp>
namespace MC
{

  template <ParticleModel _Model> class BaseParticle
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

    KOKKOS_INLINE_FUNCTION void init(Kokkos::Random_XorShift64_Pool<> _rng)
    {
      data.init(properties,_rng);
    }

    KOKKOS_INLINE_FUNCTION void
    update(double d_t, const LocalConcentrationView &concentration,MC::KPRNG _rng)
    {
      data.update(d_t, properties, concentration,_rng);
    }

    KOKKOS_INLINE_FUNCTION BaseParticle<_Model> division()
    {
      properties.status=CellStatus::IDLE;
      auto p = data.division(properties);
      // assert(p.properties.status==properties.status && properties.status==CellStatus::IDLE);
      
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