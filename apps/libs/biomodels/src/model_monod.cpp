#include <Kokkos_Atomic.hpp>
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
  constexpr double minimal_length = critcal_division_length*0.5;
  constexpr double mu_max = 10*0.77 / 3600.;
  constexpr double tau_division_proba = 1e-7;
  constexpr double ln2 = 0.69314718056;
  constexpr double tau_metabolism = 1./mu_max;
  KOKKOS_INLINE_FUNCTION double division_gamma(double lenght)
  {

    const auto length = lenght;
    if (length <= critcal_division_length)
    {
      return 0.;
    }
    constexpr double kappa = 1.;

    // return Kokkos::pow(length / critcal_division_length, kappa);
    return length/critcal_division_length;
  }

} // namespace

namespace Models
{

  KOKKOS_FUNCTION void Monod::init(MC::ParticleDataHolder &p,
                                   Kokkos::Random_XorShift64_Pool<> _rng)
  {
    this->age = 0;
    auto generator = _rng.get_state();
    this->l = generator.drand(minimal_length,critcal_division_length);
    this->mu = generator.drand(mu_max * 0.9, mu_max * 1.1);
    _rng.free_state(generator);
    _init_only_cell_lenghtening = 0.5*critcal_division_length/ln2;
    this->contrib = 0.;
  }

  KOKKOS_FUNCTION void
  Monod::update(double d_t,
                MC::ParticleDataHolder &p,
                const LocalConcentrationView &concentration,
                MC::KPRNG _rng)
  {
    const double s = concentration(0);
    if (almost_equal(s, 0.))
    {
      return;
    }
    age += d_t;
    
    const double mu_p = mu * s / (Ks + s);

    l += d_t * (YXS * mu * _init_only_cell_lenghtening);

    mu += d_t * (1.0 / tau_metabolism) * (mu_p - mu);

    contrib = mu * s / (Ks + s);

    const double proba_div =
        (1 - Kokkos::exp(-d_t / tau_division_proba)) * division_gamma(l);

    double x = _rng.double_unfiform();
    if (x < proba_div)
    {
      p.status = MC::CellStatus::CYTOKINESIS;
    }
  }

  KOKKOS_FUNCTION Monod Monod::division(MC::ParticleDataHolder &p)
  {
    age = 0;
    const double original_lenght = l;
    l = original_lenght / 2.;

    auto n = *this;

    return n;
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