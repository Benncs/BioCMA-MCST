#include <Kokkos_Printf.hpp>
#include <common/common.hpp>
#include <mc/particles/data_holder.hpp>
#include <mc/particles/particle_model.hpp>
#include <models/model_monod.hpp>

namespace
{
  constexpr double Ks = 0.01;
  constexpr double YXS = 0.5;
  constexpr double critcal_division_length = 11e-6;
  constexpr double maximal_length = 18e-6;
  constexpr double minimal_length = 7e-6;
  const double mu_max= 0.77 / 3600.;
  KOKKOS_INLINE_FUNCTION double division_gamma(double lenght)
  {

    const auto length = lenght;
    if (length <= critcal_division_length)
    {
      return 0.;
    }
    constexpr double kappa = 1.;

    return Kokkos::pow(length / critcal_division_length, kappa);
  }

} // namespace

namespace Models
{

  KOKKOS_FUNCTION void Monod::init(MC::ParticleDataHolder &p,Kokkos::Random_XorShift64_Pool<> _rng)
  {
    this->age = 0;
    auto generator = _rng.get_state();
    this->l = critcal_division_length;//generator.drand(minimal_length,critcal_division_length); 
    this->mu =  generator.drand(mu_max*0.9, mu_max*1.1); 
    _rng.free_state(generator);

    this->contrib = 0.;
  }

  KOKKOS_FUNCTION void
  Monod::update(double d_t,
                MC::ParticleDataHolder &p,
                const LocalConcentrationView &concentration,
                Kokkos::Random_XorShift64_Pool<> _rng)
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

    const double proba_div = division_gamma(l);//(1. - Kokkos::exp(-division_gamma(l) * d_t));
    auto generator = _rng.get_state();
    double x = generator.drand(0., 1.);
    _rng.free_state(generator);
    
    if (x < proba_div)
    {
      p.status = MC::CellStatus::CYTOKINESIS;
    }
  }

  KOKKOS_FUNCTION Monod Monod::division(MC::ParticleDataHolder &p)
  {
    age= 0;
    const double original_lenght = l;
    l = original_lenght / 2.;

    auto n = *this;
    
    return n;
  }

  KOKKOS_FUNCTION void Monod::contribution(MC::ParticleDataHolder &p,
                                           ContributionView contribution)
  {
    contribution(0, p.current_container) -= contrib * p.weight;
  }

  model_properties_detail_t Monod::get_properties()
  {
    return {{"mu", mu}, {"length", l},{"age",age}};
  }

  static_assert(ParticleModel<Monod>, "Check Monod Model");

} // namespace Models