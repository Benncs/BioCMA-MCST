#ifndef __BIO_MODEL_INTERDIVISION_TIME__
#define __BIO_MODEL_INTERDIVISION_TIME__

#include "mc/prng/prng.hpp"
#include <Kokkos_Core.hpp>
#include <mc/particles/particle_model.hpp>
namespace Models
{
  class InterdivisionTime
  {
    double age;
    double mu;
    double l;
    double contrib;
    double interdivision_time;

  public:
    KOKKOS_FUNCTION void init(MC::ParticleDataHolder &p,
                              MC::KPRNG _rng);

    KOKKOS_FUNCTION void update(double d_t,
                                MC::ParticleDataHolder &p,
                                const LocalConcentrationView &concentration,
                                MC::KPRNG _rng);

    KOKKOS_FUNCTION InterdivisionTime division(MC::ParticleDataHolder &p,MC::KPRNG);

    KOKKOS_FUNCTION void contribution(MC::ParticleDataHolder &p,
                                      const ContributionView& contri);

                                      KOKKOS_INLINE_FUNCTION double mass()const{return 1.;}

    model_properties_detail_t get_properties();
  };
} // namespace Models

#endif