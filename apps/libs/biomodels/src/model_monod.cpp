#include "models/utils.hpp"
#include <Kokkos_Atomic.hpp>
#include <Kokkos_Printf.hpp>
#include <common/common.hpp>
#include <mc/particles/data_holder.hpp>
#include <mc/particles/particle_model.hpp>
#include <models/model_monod.hpp>

namespace
{
  constexpr double Ks = 0.01;
  constexpr double critcal_division_length = 11e-6;
  constexpr double maximal_length = 18e-6;
  constexpr double critcal_division_length_half = critcal_division_length * 0.5;
  constexpr double minimal_length = 4.e-6;
  constexpr double mu_max = 0.77 / 3600.;
  constexpr double tau_division_proba = 1e-7;
  constexpr double ln2 = 0.69314718056;
  constexpr double tau_metabolism = (1. / mu_max);

  KOKKOS_INLINE_FUNCTION double division_gamma(double lenght)
  {
    if (lenght < critcal_division_length)
    {
      return 0.;
    }
    return lenght / critcal_division_length;
  }

} // namespace

namespace Models
{

  KOKKOS_FUNCTION void Monod::init(MC::ParticleDataHolder &p, MC::KPRNG _rng)
  {
    this->age = 0;
    auto generator = _rng.random_pool.get_state();
    // this->l = generator.drand(minimal_length, maximal_length);
    this->l = Kokkos::max(generator.normal(critcal_division_length_half,
                                           critcal_division_length_half / 5.),
                          minimal_length);
    this->mu = Kokkos::max(generator.normal(mu_max, mu_max / 5.), 0.);
    _rng.random_pool.free_state(generator);
    _init_only_cell_lenghtening =
        (critcal_division_length - critcal_division_length_half) / ln2;
    this->contrib = 0.;
  }

  KOKKOS_FUNCTION void
  Monod::update(double d_t,
                MC::ParticleDataHolder &p,
                const LocalConcentrationView &concentration,
                MC::KPRNG _rng)
  {
    double s = concentration(0);
    if (almost_equal(s, 0.))
    {
      s = 0.;
    }

    age += d_t;

    const double mu_p = mu * s / (Ks + s);

    l += d_t * (std::min(mu, mu_p) * _init_only_cell_lenghtening);

    mu += d_t * (1.0 / tau_metabolism) * (mu_p - mu);

    contrib = mu * s / (Ks + s);

    if (Models::check_probability_division(d_t, division_gamma(l), _rng))
    {
      p.status = MC::CellStatus::CYTOKINESIS;
    }
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