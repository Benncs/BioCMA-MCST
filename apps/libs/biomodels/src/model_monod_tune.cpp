#include "models/utils.hpp"
#include <Kokkos_Atomic.hpp>
#include <Kokkos_Macros.hpp>
#include <Kokkos_Printf.hpp>
#include <common/common.hpp>
#include <mc/particles/data_holder.hpp>
#include <mc/particles/particle_model.hpp>
#include <models/model_monod_tune.hpp>

constexpr double ln2 = 0.69314718056;

constexpr Models::MonodParameterSet default_set = {
    0.01, 18e-6, 11e-6, 0.5e-6, 0.77 / 3600., 1.};

namespace Models
{
  MonodParameterSet MonodTune::parameters = default_set;
  void MonodTune::set_parameters(MonodParameterSet &&values)
  {
    MonodTune::parameters = values;
  }

  KOKKOS_FUNCTION void MonodTune::init(MC::ParticleDataHolder &p,
                                       MC::KPRNG _rng)
  {
    auto generator = _rng.random_pool.get_state();

    this->l = Kokkos::max(
        MonodTune::parameters.minimal_length,
        Kokkos::max(generator.normal(MonodTune::parameters.l_0,
                                     MonodTune::parameters.l_0 / 5.),
                    0.));

    this->mu = Kokkos::max(generator.normal(MonodTune::parameters.mu_max / 4.,
                                            MonodTune::parameters.mu_max / 4),
                           0.);
    _rng.random_pool.free_state(generator);
    _init_only_cell_lenghtening =
        (MonodTune::parameters.l_1 - MonodTune::parameters.l_0) / ln2;
    this->contrib = 0.;
  }

  KOKKOS_FUNCTION void
  MonodTune::update(double d_t,
                    MC::ParticleDataHolder &p,
                    const LocalConcentrationView &concentration,
                    MC::KPRNG _rng)
  {
    const double s = Kokkos::max(0., concentration(0));
    const double mu_p =
        MonodTune::parameters.mu_max * s / (MonodTune::parameters.Ks + s);
    const double mu_eff = Kokkos::min(mu, mu_p);

    l += d_t * (mu_eff * _init_only_cell_lenghtening);

    mu += d_t * (1.0 / MonodTune::parameters.tau_metabolism) * (mu_p - mu);

    contrib = mu_eff * s / (MonodTune::parameters.Ks + s);

    Models::update_division_status(
        p.status,
        d_t,
        GammaDivision::threshold_linear(
            l, MonodTune::parameters.l_0, MonodTune::parameters.l_1),
        _rng);
  }

  KOKKOS_FUNCTION MonodTune MonodTune::division(MC::ParticleDataHolder &p)
  {
    const double original_lenght = l;

    l = original_lenght / 2.;

    return *this;
  }

  KOKKOS_FUNCTION void MonodTune::contribution(MC::ParticleDataHolder &p,
                                               ContributionView contribution)
  {
    // contribution(0, p.current_container) -= contrib * p.weight;

    auto access_contribs = contribution.access();

    // Kokkos::atomic_add(&contribution(0, p.current_container),
    //                          -contrib * p.weight);

    access_contribs(0, p.current_container) += (-contrib * p.weight);
  }

  model_properties_detail_t MonodTune::get_properties()
  {
    return {{"mu", mu}, {"lenght", l}};
  }

  static_assert(ParticleModel<MonodTune>, "Check MonodTune Model");

} // namespace Models