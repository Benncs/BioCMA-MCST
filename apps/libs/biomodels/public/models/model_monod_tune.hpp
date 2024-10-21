#ifndef __BIO_MODEL_MONOD_TUNE__
#define __BIO_MODEL_MONOD_TUNE__

#include "mc/prng/prng.hpp"
#include <Kokkos_Core.hpp>
#include <mc/particles/particle_model.hpp>
namespace Models
{
  struct MonodParameterSet
  {
    double Ks;
    double l_1;
    double l_0;
    double minimal_length;
    double mu_max;
    double tau_metabolism;
  };

  class MonodTune
  {
  private:
    // Internal properties
    double mu;
    double l;
    double contrib;
    double _init_only_cell_lenghtening;
    static MonodParameterSet parameters;

  public:
    static void set_parameters(MonodParameterSet &&values);

    KOKKOS_FUNCTION void init(MC::ParticleDataHolder &p, MC::KPRNG _rng);

    KOKKOS_FUNCTION void update(double d_t,
                                MC::ParticleDataHolder &p,
                                const LocalConcentrationView &concentration,
                                MC::KPRNG _rng);

    KOKKOS_FUNCTION MonodTune division(MC::ParticleDataHolder &p,MC::KPRNG);

    KOKKOS_FUNCTION void contribution(MC::ParticleDataHolder &p,
                                      ContributionView contri);
                                      KOKKOS_INLINE_FUNCTION double mass()const{return 1.;}

    model_properties_detail_t get_properties();
  };
} // namespace Models

#endif //__BIO_MODEL_MONOD_TUNE__