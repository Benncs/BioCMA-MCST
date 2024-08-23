#ifndef __MC_PARTICLESHPP__
#define __MC_PARTICLESHPP__

#include "common/execinfo.hpp"
#include "mc/prng/prng.hpp"
#include <cstddef>
#include <cstdint>
#include <mc/particles/data_holder.hpp>
namespace MC
{

  template <typename T>
  concept ParticleModel = requires(T model, ParticleDataHolder &p) {
    { model.init(p) } -> std::same_as<void>;
    { model.update(p) } -> std::same_as<void>;
    { model.division(p) } -> std::same_as<void>;
    { model.contribution(p) } -> std::same_as<void>;
  };

  template <ParticleModel _Model>
  class alignas(ExecInfo::cache_line_size) BaseParticle
  {
  public:

    using Model = _Model;

    KOKKOS_INLINE_FUNCTION explicit BaseParticle(double _weight = 0) noexcept
        : properties(_weight){};

    void clearState(MC::CellStatus _status = CellStatus::IDLE) noexcept
    {
      properties.reset();
      properties.status = _status;
    }

    DEFAULT_COPY_MOVE_AC(BaseParticle<_Model>)
    ~BaseParticle() = default;

    ParticleDataHolder properties;
    _Model data;
  };

  template <ParticleModel Model> using Particle = BaseParticle<Model>;

} // namespace MC

#endif