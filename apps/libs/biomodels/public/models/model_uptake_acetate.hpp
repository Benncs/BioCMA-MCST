#ifndef __BIO_MODEL_UPTAKE_ACETATE_HPP__
#define __BIO_MODEL_UPTAKE_ACETATE_HPP__

#include "mc/prng/prng.hpp"
#include <Kokkos_Core.hpp>
#include <common/kokkos_vector.hpp>
#include <mc/particles/particle_model.hpp>

namespace Models
{

  struct UptakeAcetate
  {
    using PhiUptakes = struct Contribs
    {
      double glucose;
      double acetate;
      double oxygen;
    };

    struct Rates
    {
      double glucose;
      double acetate;
      double oxygen;
      double carbon_dioxide;
      double nu;
    };

    static constexpr std::size_t n_meta = 5;

    enum class Indices : std::size_t
    {
      NU = 4,
      GLUCOSE = 0,
      O2 = 1,
      Ac = 2,
      CO2 = 3
    };

    double lenght;
    double nu;
    double a_pts;
    double a_permease;
    double n_permease;
    Rates rates;
    PhiUptakes phi_uptakes;

    KOKKOS_FUNCTION void init(MC::ParticleDataHolder &p, MC::KPRNG _rng);

    KOKKOS_FUNCTION void metabolism();

    KOKKOS_FUNCTION double uptake(const LocalConcentrationView &concentration);

    KOKKOS_FUNCTION void update(double d_t, MC::ParticleDataHolder &p, const LocalConcentrationView &concentration, MC::KPRNG _rng);

    KOKKOS_FUNCTION UptakeAcetate division(MC::ParticleDataHolder &p, MC::KPRNG _rng);

    KOKKOS_FUNCTION void contribution(MC::ParticleDataHolder &p, ContributionView contribution);

    [[nodiscard]] KOKKOS_FUNCTION double mass() const noexcept;

    model_properties_detail_t get_properties();

    template <class Archive> void serialize(Archive &ar)
    {
      ar(lenght, nu, a_pts, a_permease);
    }
  };
} // namespace Models

#endif