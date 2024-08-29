#ifndef __BIO_MODEL_MONOD__
#define __BIO_MODEL_MONOD__

#include <Kokkos_Core.hpp>
#include <mc/particles/particle_model.hpp>
namespace Models
{
  class Monod
  {
    double age;
    double mu;
    double l;
    double contrib;

  public:
    KOKKOS_FUNCTION void init(MC::ParticleDataHolder &p,Kokkos::Random_XorShift64_Pool<> _rng);

    KOKKOS_FUNCTION void update(double d_t,
                                MC::ParticleDataHolder &p,
                                const LocalConcentrationView &concentration,
                                Kokkos::Random_XorShift64_Pool<> _rng);

    KOKKOS_FUNCTION Monod division(MC::ParticleDataHolder &p);

    KOKKOS_FUNCTION void contribution(MC::ParticleDataHolder &p,
                                      ContributionView contri);

    model_properties_detail_t get_properties();
  };
} // namespace Models

#endif