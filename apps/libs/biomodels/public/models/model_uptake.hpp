#ifndef __BIO_model_test_HPP__
#define __BIO_model_test_HPP__

#include "mc/prng/prng.hpp"
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
    std::array<double,2> contribs{0,0};
    KOKKOS_FUNCTION void init(MC::ParticleDataHolder &p, MC::KPRNG _rng);

    KOKKOS_FUNCTION void update(double d_t,
                                MC::ParticleDataHolder &p,
                                const LocalConcentrationView &concentration,
                                MC::KPRNG _rng);

    KOKKOS_FUNCTION Uptake division(MC::ParticleDataHolder &p,MC::KPRNG);

    KOKKOS_FUNCTION void contribution(MC::ParticleDataHolder &p,
                                      const ContributionView& contri);

    [[nodiscard]]  KOKKOS_FUNCTION  double mass() const noexcept;

    model_properties_detail_t get_properties();
  };
} // namespace Models

#endif //__BIO_model_test_HPP__
