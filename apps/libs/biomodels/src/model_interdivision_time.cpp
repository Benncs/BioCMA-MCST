
// #include <common/common.hpp>
#include <mc/particles/data_holder.hpp>
#include <mc/particles/particle_model.hpp>
#include <models/model_interdivision_time.hpp>

namespace
{
  constexpr double Ks = 0.01;
  constexpr double YXS = 0.5;
  constexpr double critical_division_length = 11e-6;
  constexpr double minimal_length = 7e-6;
  constexpr double mu_max = 0.77 / 3600.;
  constexpr double ln2 = 0.69314718056;
  constexpr double tau_division = ln2 / mu_max;
} // namespace

namespace Models
{

  KOKKOS_FUNCTION void InterdivisionTime::init(MC::ParticleDataHolder &p,
                                               MC::KPRNG _rng)
  {
    this->age = 0;
    this->interdivision_time = tau_division;
    this->contrib = 0.;
    {
      auto generator = _rng.random_pool.get_state();
      this->mu = generator.drand(mu_max * 0.9, mu_max * 1.1);
      this->l = generator.drand(minimal_length, critical_division_length);
      _rng.random_pool.free_state(generator);
    }
  }

  KOKKOS_FUNCTION void
  InterdivisionTime::update(double d_t,
                            MC::ParticleDataHolder &p,
                            const LocalConcentrationView &concentration,
                            MC::KPRNG /*_rng*/)
  {
    const double s = concentration(0);
    if (almost_equal(s, 0.))
    {
      return;
    }
    age += d_t;
    constexpr double cell_lenghtening = 6.51e-6;
    constexpr auto tau_metabolism = 3600;
    const double mu_p = mu * s / (Ks + s);

    l += d_t * (YXS * mu * cell_lenghtening);

    mu += d_t * (1.0 / tau_metabolism) * (mu_p - mu);

    contrib = mu * s / (Ks + s);

    if (interdivision_time < age && !almost_equal(s, 0., 1e-4))
    {
      p.status = MC::CellStatus::CYTOKINESIS;
    }
  }

  KOKKOS_FUNCTION InterdivisionTime
  InterdivisionTime::division(MC::ParticleDataHolder &p,MC::KPRNG)
  {
    age = 0;
    const double original_lenght = l;
    l = original_lenght / 2.;

    auto n = *this;

    return n;
  }

  KOKKOS_FUNCTION void
  InterdivisionTime::contribution(MC::ParticleDataHolder &p,
                                  const ContributionView& contribution)
  {
    // contribution(0, p.current_container) -= contrib * p.weight;

    auto access_contribs = contribution.access();

    // Kokkos::atomic_add(&contribution(0, p.current_container),
    //                          -contrib * p.weight);

    access_contribs(0,p.current_container)+= (-contrib*p.weight);
  }

  model_properties_detail_t InterdivisionTime::get_properties()
  {
    return {{"mu", mu},
            {"lenght", l},
            {"age", age},
            {"interdivision_time", interdivision_time}};
  }

  static_assert(ParticleModel<InterdivisionTime>,
                "Check InterdivisionTime Model");

} // namespace Models