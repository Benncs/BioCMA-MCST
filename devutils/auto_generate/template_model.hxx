#ifndef __BIO_@__model__name__@_HPP__
#  define __BIO_@__model__name__@_HPP__

#  include <mc/prng/prng.hpp>
#  include <Kokkos_Core.hpp>
#  include <common/kokkos_vector.hpp>
#  include <mc/particles/particle_model.hpp>

namespace Models
{
  struct @__model__name__@
  {
    KOKKOS_FUNCTION void init(MC::ParticleDataHolder &p, MC::KPRNG _rng);

  KOKKOS_FUNCTION void update(double d_t,
                              MC::ParticleDataHolder &p,
                              const LocalConcentrationView &concentration,
                              MC::KPRNG _rng);

  KOKKOS_FUNCTION @__model__name__@ division(MC::ParticleDataHolder &p);

  KOKKOS_FUNCTION void contribution(MC::ParticleDataHolder &p,
                                    ContributionView contri);

  model_properties_detail_t get_properties();
}; 
} // namespace Models

#endif //__BIO_@__model__name__@_HPP__

