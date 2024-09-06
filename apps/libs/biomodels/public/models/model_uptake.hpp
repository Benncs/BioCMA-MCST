#ifndef __BIO_model_test_HPP__
#define __BIO_model_test_HPP__

#include "mc/prng/prng.hpp"
#include <Kokkos_Core.hpp>
#include <common/kokkos_vector.hpp>
#include <mc/particles/particle_model.hpp>

namespace Models
{
  struct Uptake
  {
    double lenght;
    double nu;
    double a_pts;
    double a_permease;
    double n_permease;
    double _init_only_cell_lenghtening;
    double contrib;
    KOKKOS_FUNCTION void init(MC::ParticleDataHolder &p, MC::KPRNG _rng);

    KOKKOS_FUNCTION void update(double d_t,
                                MC::ParticleDataHolder &p,
                                const LocalConcentrationView &concentration,
                                MC::KPRNG _rng);

    KOKKOS_FUNCTION Uptake division(MC::ParticleDataHolder &p);

    KOKKOS_FUNCTION void contribution(MC::ParticleDataHolder &p,
                                      ContributionView contri);

    model_properties_detail_t get_properties();
  };
} // namespace Models

#endif //__BIO_model_test_HPP__
