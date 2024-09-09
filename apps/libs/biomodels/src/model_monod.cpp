#include "models/utils.hpp"
#include <Kokkos_Atomic.hpp>
#include <Kokkos_Macros.hpp>
#include <Kokkos_Printf.hpp>
#include <common/common.hpp>
#include <mc/particles/data_holder.hpp>
#include <mc/particles/particle_model.hpp>
#include <models/model_monod.hpp>

namespace
{
  constexpr double Ks = 0.01;
  constexpr double l_1 = 18e-6;
  constexpr double l_0 = 11e-6;
  constexpr double minimal_length = 0.5e-6;
  constexpr double mu_max = 0.77 / 3600.;
  constexpr double ln2 = 0.69314718056;
  constexpr double tau_metabolism = (1. / mu_max);

  KOKKOS_INLINE_FUNCTION double division_gamma(double lenght)
  {
    static constexpr double denum = l_1 - l_0;
    return (lenght < l_0) ? 0 : (lenght - l_0) / denum;
  }

} // namespace

namespace Models
{

  KOKKOS_FUNCTION void Monod::init(MC::ParticleDataHolder &p, MC::KPRNG _rng)
  {
    this->age = 0;
    auto generator = _rng.random_pool.get_state();

    this->l = Kokkos::max(minimal_length,
                          Kokkos::max(generator.normal(l_0, l_0 / 5.), 0.));

    this->mu = Kokkos::max(generator.normal(mu_max/2., mu_max / 5), 0.);
    _rng.random_pool.free_state(generator);
    _init_only_cell_lenghtening = (l_1 - l_0) / ln2;
    this->contrib = 0.;
  }

  KOKKOS_FUNCTION void
  Monod::update(double d_t,
                MC::ParticleDataHolder &p,
                const LocalConcentrationView &concentration,
                MC::KPRNG _rng)
  {
    const double s = Kokkos::max(0., concentration(0));
    const double mu_p = mu_max * s / (Ks + s);
    const double mu_eff = Kokkos::min(mu, mu_p);

    age += d_t;

    l += d_t * (mu_eff * _init_only_cell_lenghtening);

    mu += d_t * (1.0 / tau_metabolism) * (mu_p - mu);

    contrib = mu_eff * s / (Ks + s);

    Models::update_division_status(
        p.status, d_t, GammaDivision::threshold_linear(l, l_0, l_1), _rng);
  }

  KOKKOS_FUNCTION Monod Monod::division(MC::ParticleDataHolder &p)
  {
    age = 0;
    const double original_lenght = l;

    l = original_lenght / 2.;

    return *this;
  }

  KOKKOS_FUNCTION void Monod::contribution(MC::ParticleDataHolder &p,
                                           ContributionView contribution)
  {
    // contribution(0, p.current_container) -= contrib * p.weight;

    Kokkos::atomic_add_fetch(&contribution(0, p.current_container),
                             -contrib * p.weight);
  }

  model_properties_detail_t Monod::get_properties()
  {
    return {{"mu", mu}, {"lenght", l}, {"age", age}};
  }

  static_assert(ParticleModel<Monod>, "Check Monod Model");

} // namespace Models