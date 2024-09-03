#include "common/common.hpp"
#include "mc/prng/prng.hpp"
#include "models/utils.hpp"
#include <Kokkos_Macros.hpp>
#include <models/model_uptake.hpp>

static constexpr double tauPTS = 25.0;
static constexpr double tau_f = 5.0;
static constexpr double tau_d = 5.0;
static constexpr double NPermease_max = 5e4;
static constexpr double NPermease_init = 1e4;
static constexpr double tauAu = 5.0;
static constexpr double tauAd = 5.0;
static constexpr double tau_metabolism = 3600;
static constexpr double phi_pts_max = 4.454e-12 / 3600.;
static constexpr double kpts = 1e-2;
static constexpr double kppermease = 1e-3;
static constexpr double psi_permease = 1.25e-13 / 3600.;
static constexpr double YXS = 0.5;
static constexpr double YXO = 1e-4;
static constexpr double YSA = 1e-4;
static constexpr double psi_o_meta = 20e-3 / 3600 * 32e-3; // 20mmmol O2 /h;
static constexpr double critcal_division_length = 11e-6;
static constexpr double critcal_division_length_half =
    0.5 * critcal_division_length;
static constexpr double minimal_length = 7e-6;
static constexpr double maximal_length = 18e-6;
static constexpr double mu_max = 10.*0.77 / 3600.;
static constexpr double ln2 = 0.69314718056;

KOKKOS_INLINE_FUNCTION double phi_pts(double a_pts, double S)
{
  return a_pts * (phi_pts_max * S / (kpts + S));
}

KOKKOS_INLINE_FUNCTION double
phi_permease(double n_permease, double a_permease, double S)
{

  return n_permease * psi_permease * a_permease *
         (phi_pts_max * S / (kppermease + S));
}

namespace Models
{
  KOKKOS_FUNCTION void Uptake::init(MC::ParticleDataHolder &p, MC::KPRNG _rng)
  {
    auto generator = _rng.random_pool.get_state();
    // this->l = generator.drand(minimal_length, maximal_length);
    this->lenght =
        Kokkos::max(generator.normal(critcal_division_length_half,
                                     critcal_division_length_half / 3.),
                    minimal_length);

    this->a_permease = 1e-3;
    this->a_pts = generator.drand(0., 1.);
    this->n_permease = generator.drand(0., 1.) * NPermease_init;
    this->mu_eff = Kokkos::max(generator.normal(mu_max, mu_max / 5.), 0.);
    _rng.random_pool.free_state(generator);

    _init_only_cell_lenghtening =
        (critcal_division_length - critcal_division_length_half) / ln2;
    this->contrib = 0.;
  }

  KOKKOS_FUNCTION void
  Uptake::update(double d_t,
                 MC::ParticleDataHolder &p,
                 const LocalConcentrationView &concentration,
                 MC::KPRNG _rng)
  {
    double s = concentration(0);
    if (almost_equal(s, 0.))
    {
      s = 0.;
    }
    const double phi_s_pts = phi_pts(a_pts, s);
    const double phi_s_in = phi_s_pts + phi_permease(n_permease, a_permease, s);
    const double gamma_PTS_S = phi_s_pts / phi_pts_max;
    const double mu_p = YXS * phi_s_in;
  
    this->lenght += d_t * (mu_eff * _init_only_cell_lenghtening);
    this->mu_eff += d_t * ((1.0 / tau_metabolism) * (mu_p - this->mu_eff));
    this->a_pts += d_t * ((1.0 / tauPTS) * ((s / (kpts + s)) - this->a_pts));

    this->a_permease +=
        d_t *
        (((1.0 / tauAu) * gamma_PTS_S + (1.0 / tauAd) * (1.0 - gamma_PTS_S)) *
         (1.0 - gamma_PTS_S - a_permease));

    this->n_permease +=
        d_t * ((1.0 / tau_f) * (kppermease / (kppermease + s)) +
               (1.0 / tau_d) * (s / (kpts + s)) *
                   (NPermease_max * (1.0 - gamma_PTS_S) - a_permease));

    this->contrib = phi_pts(a_pts, s) + phi_permease(n_permease, a_permease, s);

    const double gamma = (lenght > critcal_division_length) ? 1 : 0;
    if (Models::check_probability_division(d_t, gamma, _rng))
    {
      p.status = MC::CellStatus::CYTOKINESIS;
    }
  }

  KOKKOS_FUNCTION Uptake Uptake::division(MC::ParticleDataHolder &p)
  {
    const double original_lenght = this->lenght;
    const double original_n_permease = this->n_permease;

    this->lenght = original_lenght / 2.;
    this->n_permease = original_n_permease / 2;
    return *this;
  }

  KOKKOS_FUNCTION void Uptake::contribution(MC::ParticleDataHolder &p,
                                            ContributionView contribution)
  {

    Kokkos::atomic_add_fetch(&contribution(0, p.current_container),
                             -this->contrib * p.weight);
  }

  model_properties_detail_t Uptake::get_properties()
  {
    return {{"lenght", lenght},
            {"a_pts", a_pts},
            {"mu_eff", mu_eff},
            {"a_permease", a_permease},
            {"n_permease", n_permease}};
  }

}; // namespace Models
