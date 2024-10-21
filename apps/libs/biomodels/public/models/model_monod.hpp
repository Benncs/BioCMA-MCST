#ifndef __BIO_MODEL_MONOD__
#define __BIO_MODEL_MONOD__

#include "common/kokkos_vector.hpp"
#include "mc/prng/prng.hpp"
#include <Kokkos_Core.hpp>
#include <mc/particles/particle_model.hpp>
namespace Models
{

  class alignas(ExecInfo::cache_line_size) Monod
  {
    double mu;
    double l;
    double contrib;

  private:
    double _init_only_cell_lenghtening;

  public:
    KOKKOS_FUNCTION void init(MC::ParticleDataHolder &p, MC::KPRNG _rng)noexcept;

    KOKKOS_FUNCTION void update(double d_t,
                                MC::ParticleDataHolder &p,
                                const LocalConcentrationView &concentration,
                                MC::KPRNG _rng)noexcept;

    KOKKOS_FUNCTION Monod division(MC::ParticleDataHolder &p,MC::KPRNG)noexcept;

    KOKKOS_FUNCTION void contribution(MC::ParticleDataHolder &p,
                                      ContributionView contri)noexcept;
    KOKKOS_FUNCTION [[nodiscard]] double mass() const noexcept;

    model_properties_detail_t get_properties() noexcept;
  };
} // namespace Models

#endif