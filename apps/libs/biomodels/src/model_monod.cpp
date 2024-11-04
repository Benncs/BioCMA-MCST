#include <Kokkos_Atomic.hpp>
#include <Kokkos_Macros.hpp>
#include <Kokkos_Printf.hpp>
#include <common/common.hpp>
#include <mc/particles/data_holder.hpp>
#include <mc/particles/particle_model.hpp>
#include <models/model_monod.hpp>
#include <models/utils.hpp>

namespace
{
  constexpr double y_s_x = 2.;
  constexpr double Ks = 0.01;
  constexpr double l_1 = 1.7e-6;
  constexpr double l_0 = 0.9e-6;
  constexpr double minimal_length = 0.5e-6;
  constexpr double mu_max = 0.77 / 3600.;
  constexpr double ln2 = 0.69314718056;
  constexpr double tau_metabolism = (1. / mu_max);
} // namespace
static constexpr double factor =  1000 * 3.14 * (0.8e-6) * (0.8e-6) / 4.;
#define MASS() (factor * l)

namespace Models
{

  KOKKOS_FUNCTION void Monod::init(MC::ParticleDataHolder &p, MC::KPRNG _rng) noexcept
  {
    auto generator = _rng.random_pool.get_state();

    this->l = l_0 / 2.;
    // Kokkos::max(minimal_length,
    //             Kokkos::max(generator.normal(l_0 / 2., l_0 / 2. / 5.), 0.));

    this->mu = Kokkos::max(generator.normal(mu_max/2. , mu_max / 4), 0.);
    _rng.random_pool.free_state(generator);
    static_assert(l_1 > l_0, "Monod: Bad Model Parameter ");
    _init_only_cell_lenghtening =l_0/2. / ln2;

    // p.weight = p.weight/mass();
    this->contrib = 0.;
  }

  KOKKOS_FUNCTION void
  Monod::update(double d_t,
                MC::ParticleDataHolder &p,
                const LocalConcentrationView &concentration,
                MC::KPRNG _rng) noexcept
  {
    const double s = Kokkos::max(0., concentration(0));
    const double mu_p = mu_max * s / (Ks + s);
    const double mu_eff = Kokkos::min(mu, mu_p);
    contrib = mu_eff *y_s_x*MASS() ;
    l += d_t * (mu_eff * _init_only_cell_lenghtening);
    mu += d_t * (1.0 / tau_metabolism) * (mu_p - mu);

    Models::update_division_status(
        p.status, d_t, GammaDivision::threshold_linear(l, l_0, l_1), _rng);
  }

  KOKKOS_FUNCTION Monod Monod::division(MC::ParticleDataHolder &p,MC::KPRNG) noexcept
  {
    const double original_lenght = l;

    l = original_lenght / 2.;

    return *this;
  }

  KOKKOS_FUNCTION void Monod::contribution(MC::ParticleDataHolder &p,
                                           ContributionView contribution) noexcept
  {
    // contribution(0, p.current_container) -= contrib * p.weight;
    auto access_contribs = contribution.access();

    // Kokkos::atomic_add(&contribution(0, p.current_container),
    //                          -contrib * p.weight);

    access_contribs(0, p.current_container) += (-contrib * p.weight);

    // Kokkos::atomic_add_fetch(&contribution(0, p.current_container),
    //                          -contrib * p.weight);

    //                          Kokkos::atomic_add(&contribution(0,
    //                          p.current_container), -contrib * p.weight);
  }

  model_properties_detail_t Monod::get_properties() noexcept
  {
    return {{"mu", mu}, {"lenght", l}, {"mass", mass()}};
  }

  KOKKOS_FUNCTION double Monod::mass() const noexcept
  {

    return MASS();
  }

  static_assert(ParticleModel<Monod>, "Check Monod Model");

} // namespace Models
