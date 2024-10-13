#include "common/common.hpp"
#include "mc/prng/prng.hpp"
#include "models/utils.hpp"
#include <Kokkos_Macros.hpp>
#include <models/model_uptake.hpp>

static constexpr double tauPTS = 50.0;
static constexpr double tau_f = 100;
static constexpr double tau_d = 6 * 60;
static constexpr double NPermease_max = 200;
static constexpr double NPermease_init = NPermease_max;
static constexpr double tauAu = 5.;
static constexpr double tauAd = 15.;

static constexpr double phi_pts_max = 5.1451e-12/ 3600. * 1e-3;
static constexpr double ks = 1e-2;
static constexpr double kpts = 1e-3;
static constexpr double kppermease = 1e-2;
static constexpr double psi_permease = 1.28625e-13 / 3600. * 1e-3; // par permease
static constexpr double YXS = 2.;
static constexpr double YXO = 1e-4;
static constexpr double YSA = 1e-4;
static constexpr double psi_o_meta = 20e-3 / 3600 * 32e-3; // 20mmmol O2 /h;
static constexpr double mu_max = 0.77 / 3600.;
static constexpr double tau_metabolism = 1. / mu_max;
static constexpr double ln2 = 0.69314718056;
static constexpr double l_1 = 1.7e-6;
static constexpr double l_0 = 0.9e-6;
constexpr double minimal_length = 0.4e-6;

constexpr double y = phi_pts_max/(psi_permease*NPermease_max);


#define MONOD_RATIO(__c1__, __x1__, __k1__) ((__c1__) * (__x1__) / ((__x1__) + (__k1__)))

KOKKOS_INLINE_FUNCTION double phi_pts(double a_pts, double S)
{
  return a_pts * MONOD_RATIO(phi_pts_max, S, kpts);
}

KOKKOS_INLINE_FUNCTION double
phi_permease(double n_permease, double a_permease, double S)
{
  return a_permease * MONOD_RATIO(n_permease * psi_permease, S, kppermease);
}

namespace Models
{
  KOKKOS_FUNCTION void Uptake::init(MC::ParticleDataHolder &p, MC::KPRNG _rng)
  {
    auto generator = _rng.random_pool.get_state();

    // this->lenght = Kokkos::max(
    //     minimal_length, Kokkos::max(generator.normal(l_0/2, l_0 / 5.), 0.));
    this->lenght = l_0 / 2.;
    this->a_permease = Kokkos::max(generator.normal(1e-3, 1e-4), 0.);

    this->a_pts = Kokkos::min(1., Kokkos::max(generator.normal(0.15, 0.1), 0.));
    this->n_permease = Kokkos::max(generator.normal(NPermease_init / 2., NPermease_init / 5.), 0.);
    
    this->nu = Kokkos::max(generator.normal(1e-6/3600, 1e-6/3600/5. / 5.), 0.);
    _rng.random_pool.free_state(generator);

    _init_only_cell_lenghtening = l_0 / 2. / ln2;
    this->contrib = 0.;
  }

  KOKKOS_FUNCTION void
  Uptake::update(double d_t,
                 MC::ParticleDataHolder &p,
                 const LocalConcentrationView &concentration,
                 MC::KPRNG _rng)
  {

       //dmdt = phis_s_in/YXS
    //pid2/4prhodl*dl/dt=phi_s_in/YXS
    //dl/dt = phis_in/yxs/

    static constexpr double factor = 1000.*3.14*(0.8e-6)*(0.8e-6)/4.;

    const double s = Kokkos::max(0., concentration(0));
    const double phi_s_pts = phi_pts(a_pts, s);
    const double phi_s_in = phi_s_pts + phi_permease(n_permease, a_permease, s);
    const double gamma_PTS_S = phi_s_pts / phi_pts_max;

    const double nu_eff = Kokkos::min(phi_s_in / factor / YXS, nu);

   
    this->lenght += d_t * nu_eff;
    this->nu += d_t * (1.0 / tau_metabolism) * (nu_eff - this->nu);
    this->a_pts += d_t * 1.0 / tauPTS * (MONOD_RATIO(1., s, kpts) - this->a_pts);

    this->a_permease +=
        d_t *
        (((1.0 / tauAu) * gamma_PTS_S + (1.0 / tauAd) * (1.0 - gamma_PTS_S)) *
         (1.0 - gamma_PTS_S - a_permease));

    this->n_permease += d_t *
                        (MONOD_RATIO(1. / tau_f, kppermease, s) + MONOD_RATIO(1. / tau_d, s, kppermease)) * (MONOD_RATIO(NPermease_max, ks, s) - n_permease);

    this->contrib = phi_s_in;

    Models::update_division_status(
        p.status, d_t, GammaDivision::threshold_linear(lenght, l_0, l_1), _rng);
  }

  KOKKOS_FUNCTION Uptake Uptake::division(MC::ParticleDataHolder &p)
  {
    const double original_lenght = this->lenght;
    const double original_n_permease = this->n_permease;

    this->lenght = original_lenght / 2.;
    this->n_permease = original_n_permease / 2;
    return *this;
  }

  KOKKOS_FUNCTION double Uptake::mass() const noexcept
  {
    return (1000 * 3.14 * (0.8e-6) * (0.8e-6) / 4. * lenght);
  }

  KOKKOS_FUNCTION void Uptake::contribution(MC::ParticleDataHolder &p,
                                            ContributionView contribution)
  {

    auto access_contribs = contribution.access();

    // Kokkos::atomic_add(&contribution(0, p.current_container),
    //                          -contrib * p.weight);

    access_contribs(0, p.current_container) += -contrib * p.weight;
  }

  model_properties_detail_t Uptake::get_properties()
  {
    return {{"lenght", lenght},
            {"a_pts", a_pts},
            {"mu_eff", nu},
            {"a_permease", a_permease},
            {"n_permease", n_permease},
            {"mass", mass()}};
  }

}; // namespace Models
