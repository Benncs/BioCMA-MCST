#ifndef __BIO_MODEL_MONOD__
#define __BIO_MODEL_MONOD__

#include "Kokkos_Assert.hpp"
#include "Kokkos_Macros.hpp"
#include "Kokkos_Printf.hpp"
#include "Kokkos_UnorderedMap.hpp"
#include "common/kokkos_vector.hpp"
#include "mc/prng/prng.hpp"
#include <Kokkos_Core.hpp>
#include <Kokkos_Random.hpp>
#include <cstddef>
#include <mc/particles/particle_model.hpp>
#include <models/utils.hpp>
#include <string>
namespace implMonod
{
  constexpr double y_s_x = 2.;
  constexpr double Ks = 0.01;
  constexpr double l_1 = 5e-6;
  constexpr double l_0 = 2e-6;
  constexpr double minimal_length = 0.5e-6;
  constexpr double mu_max = 0.77 / 3600.;
  constexpr double ln2 = 0.69314718056;
  constexpr double tau_metabolism = (1. / mu_max);
} // namespace implMonod

static constexpr double factor = 1000 * 3.14 * (0.8e-6) * (0.8e-6) / 4.;
#define MASS() (factor * l)
namespace Models
{

  class Monod
  {
    double mu;
    double l;
    double contrib;

  private:
    double _init_only_cell_lenghtening;

  public:

  
    KOKKOS_FUNCTION void init(MC::ParticleDataHolder& p, MC::KPRNG _rng) noexcept
    {
      using namespace implMonod;
      constexpr double local_l = minimal_length;
      auto generator = _rng.random_pool.get_state();

      // this->l = 4e-6 *0.9;//l_0 / 2.;
      this->l = Kokkos::max(local_l, Kokkos::max(generator.normal(l_0 * 2., l_0 * 2. / 5.), 0.));

      this->mu = Kokkos::max(generator.normal(mu_max / 2., mu_max / 4), 0.);
      _rng.random_pool.free_state(generator);
      static_assert(l_1 > l_0, "Monod: Bad Model Parameter ");
      _init_only_cell_lenghtening = l_0 / 2. / ln2;

      // p.weight = p.weight/mass();
      this->contrib = 0.;
    }

    KOKKOS_FUNCTION void update(double d_t,
                                MC::ParticleDataHolder& p,
                                const LocalConcentrationView& concentration,
                                Kokkos::Random_XorShift64_Pool<> _rng) noexcept
    {
      using namespace implMonod;
      const double s = Kokkos::max(0., concentration(0));
      const double mu_p = mu_max * s / (Ks + s);
      const double mu_eff = Kokkos::min(mu, mu_p);
      contrib = mu_eff * y_s_x * MASS();
      l += d_t * (mu_eff * _init_only_cell_lenghtening);
      mu += d_t * (1.0 / tau_metabolism) * (mu_p - mu);
      const auto gamma = GammaDivision::threshold_linear(l, l_0, l_1);
      Models::update_division_status(p.status, d_t, gamma, _rng);
    }

    KOKKOS_FUNCTION Monod division(MC::ParticleDataHolder& p, MC::KPRNG) noexcept
    {
      using namespace implMonod;
      const double original_lenght = l;

      l = original_lenght / 2.;
      auto child = *this;
      KOKKOS_ASSERT(child.l==this->l && child.l!=0.);
      return child;
    }

    KOKKOS_FUNCTION void contribution(MC::ParticleDataHolder& p,
                                      const ContributionView& contribution) noexcept
    {
      using namespace implMonod;
      auto access_contribs = contribution.access();

      access_contribs(0, p.current_container) += (-contrib * p.weight);
    }
    
    [[nodiscard]] KOKKOS_FUNCTION double mass() const noexcept
    {
      return MASS();
    }

    template <class Archive> void serialize(Archive& ar)
    {
      using namespace implMonod;
      ar(mu, l, contrib, _init_only_cell_lenghtening);
    }

    // model_properties_detail_t get_properties() noexcept
    // {
    //   return {{"mu", mu}, {"length", l}, {"mass", mass()}};
    // }

    static std::vector<std::string> names()
    {
      using namespace implMonod;
      return {"mu", "length", "mass"};
    }

    KOKKOS_INLINE_FUNCTION static consteval std::size_t get_number()
    {
      return 3;
    }

    KOKKOS_INLINE_FUNCTION void fill_properties(SubViewtype full) const
    {
      full(0) = mu;
      full(1) = l;
      full(2) = MASS();
    }
  };
} // namespace Models

#endif
