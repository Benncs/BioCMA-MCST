#ifndef __BIO_MODEL_UPTAKE_ACETATE_HPP__
#define __BIO_MODEL_UPTAKE_ACETATE_HPP__

#include "mc/prng/prng.hpp"
#include "models/utils.hpp"
#include <Kokkos_Core.hpp>
#include <Kokkos_Core_fwd.hpp>
#include <Kokkos_Macros.hpp>
#include <algorithm>
#include <common/kokkos_vector.hpp>
#include <mc/particles/particle_model.hpp>
namespace
{
#define MONOD_RATIO(__c1__, __x1__, __k1__) ((__c1__) * (__x1__) / ((__x1__) + (__k1__)))
#define _INDICES(__name__) static_cast<std::size_t>(Indices::__name__)

  struct Yields
  {
    static constexpr double x_s = 2.;
    static constexpr double Y_EG = 12.05 / 100; // mol_E/mol_G
    static constexpr double Y_XG = 0.8;         // mol_X/mol_G
    static constexpr double m = 250;            // µmol_G/g_X.h
    static constexpr double Y_OG = 6. / 100.;   // mol_O/mol_G
    static constexpr double Yo_EG = 20. / 100.; // mol_E/mol_G
    static constexpr double Yf_EG = 3. / 100.;  // mol_E/mol_G
    static constexpr double Y_AG = 1.5;         // mol_A/mol_G (modified for CO2 production)
    static constexpr double Y_CO2G_over = 2;    // mol_CO2/mol_G (overflow metabolism)

    static constexpr double Y_EA = 12.05 / (3 * 100); // mol_E/mol_A
    static constexpr double Y_XA = 0.4 / 100;         // mol_X/mol_A
    static constexpr double Y_OA = 2 / 100;           // mol_O/mol_A
    static constexpr double Yo_EA = 4.67 / 100;       // mol_E/mol_A
  };

  struct MolarMass
  {
    static constexpr double glucose = 180e-3;
    static constexpr double dioxygen = 32e-3;
    static constexpr double acetate = 59e-3;
    static constexpr double co2 = 44 - 3;
  };

  // All in seconds
  struct Tau
  {
    static constexpr double new_permease = 40.;
    static constexpr double tau_rm_perm = 200.;
    static constexpr double pts = 20.;
    static constexpr double Au = 40.;
    static constexpr double Ad = 5.;
  };

  constexpr double l_1 = 9e-6;
  constexpr double l_0 = 6e-6;
  constexpr double d_m = 1e-6;
  constexpr double factor = 1000. * 3.14 * d_m * d_m / 4.;
  static_assert(l_0 < l_1, "Length");

  consteval double f_num_max()
  {
    return l_0 * 0.69314718056 / 3600.; // 0.69314718056=ln(2)
  }

  consteval double f_phi_pts_max(double nu_m)
  {
    return nu_m * Yields::x_s * factor;
  }

  constexpr double nu_max = f_num_max();

  constexpr double phi_pts_max = f_phi_pts_max(nu_max);
  constexpr double phi_o2_max = 15e-6 / 3600 / 1000; // mol0/kgx/s
  constexpr double phi_a_max = phi_pts_max / 3.;     // mol0/kgx/s
  constexpr double k_pts = 1e-3;
  constexpr double kppermease = 1e-2;
  constexpr double phi_permease_specific = phi_pts_max / 40.;

  constexpr double psi_o_meta = 20e-3 / 3600 * 32e-3; // 20mmmol O2 /h;
  constexpr double tau_metabolism = 2 * 3600;

  constexpr double NPermease_max = 200;
  constexpr double NPermease_init = 0.;

  KOKKOS_INLINE_FUNCTION double phi_pts(double a_pts, double S)
  {
    return a_pts * MONOD_RATIO(phi_pts_max, S, k_pts);
  }

  KOKKOS_INLINE_FUNCTION double phi_permease(double n_permease, double a_permease, double S)
  {
    return a_permease * MONOD_RATIO(n_permease * phi_permease_specific, S, kppermease);
  }

  static_assert(NPermease_max * phi_permease_specific == 5. * phi_pts_max, "pts/perm");

} // namespace

namespace Models
{

  struct UptakeAcetate
  {
    static constexpr std::size_t n_meta = 5;
    enum class Indices : std::size_t
    {
      NU = 4,
      GLUCOSE = 0,
      O2 = 1,
      Ac = 2,
      CO2 = 3
    };

    double lenght;
    double nu;
    double a_pts;
    double a_permease;
    double n_permease;

    double contrib_s;
    double contrib_a;
    double contrib_c02;

    std::array<double, 3> contribs{};
    std::array<double, n_meta> q{};

    KOKKOS_INLINE_FUNCTION void init(MC::ParticleDataHolder &p, MC::KPRNG _rng)
    {
      auto generator = _rng.random_pool.get_state();

      // this->lenght = Kokkos::max(
      //     minimal_length, Kokkos::max(generator.normal(l_0/2, l_0 / 3.), 0.));
      this->lenght = l_0 * 0.9;
      this->a_permease = Kokkos::max(generator.normal(1e-3, 1e-4), 0.);

      this->a_pts = Kokkos::min(1., Kokkos::max(generator.normal(0.8, 0.1), 0.));
      this->n_permease = Kokkos::max(generator.normal(NPermease_init / 2., NPermease_init / 5.), 0.);

      this->nu = Kokkos::max(generator.normal(nu_max / 5., nu_max / 5. / 5.), 0.);
      _rng.random_pool.free_state(generator);
    }

    KOKKOS_INLINE_FUNCTION void metabolism(const double phi_s_in, const double phi_o2_in, double phi_a_in)
    {

      const double qr1max = phi_s_in;

      const double phi_s_in_mol = phi_s_in / (mass() * MolarMass::glucose);

      const double qmodeB = Kokkos::min((Yields::Y_EG + Yields::Yo_EG) / Yields::Yo_EG * qr1max,
                                        Kokkos::min(phi_s_in_mol, (Yields::Y_EG + Yields::Yo_EG) / (Yields::Y_EG * Yields::Y_OG) * phi_o2_in));
      const double qsr = phi_s_in_mol - qmodeB;
      const double qO2r = phi_o2_in - (Yields::Y_EG * Yields::Y_OG) / (Yields::Y_EG + Yields::Yo_EG) * qmodeB;
      const double vpr = (Yields::Y_EG + Yields::Yo_EG) / Yields::Yo_EG * qr1max - qmodeB;
      const double qmodeC = Kokkos::min((Yields::Y_EG + Yields::Yf_EG) / Yields::Yf_EG * vpr, qsr);
      const double qmodeD = Kokkos::min(qsr, qO2r);
      const double qmodeE = Kokkos::min((Yields::Y_EA + Yields::Yo_EA) / Yields::Yo_EA * vpr,
                                        Kokkos::min(phi_a_in, (Yields::Y_EA + Yields::Yo_EA) / (Yields::Yo_EA * Yields::Y_OA) * qO2r));

      q[_INDICES(GLUCOSE)] = -phi_s_in;
      q[_INDICES(O2)] = mass() * MolarMass::dioxygen *
                        (Yields::Y_EG * Yields::Y_OG / (Yields::Y_EG + Yields::Yo_EG) * qmodeB +
                         Yields::Y_EA * Yields::Y_OA / (Yields::Y_EA + Yields::Yo_EA) * qmodeE + qmodeD);
      q[_INDICES(Ac)] =
          mass() * (Yields::Y_AG * Yields::Y_EG / (Yields::Y_EG + Yields::Yf_EG) * qmodeC + Yields::Y_AG * qmodeD - qmodeE) * MolarMass::acetate;

      q[_INDICES(CO2)] = mass() * (q[_INDICES(O2)] / MolarMass::dioxygen) * MolarMass::co2;
      
      const double MMx = 113.1; // REMOVE
      q[_INDICES(NU)] = (Yields::Y_XG * Yields::Yo_EG / (Yields::Y_EG + Yields::Yo_EG) * qmodeB +
                         Yields::Y_XG * Yields::Yf_EG / (Yields::Y_EG + Yields::Yf_EG) * qmodeC +
                         Yields::Y_XA * Yields::Yo_EA / (Yields::Y_EA + Yields::Yo_EA) * qmodeE) *
                        MMx * lenght;
    }

    KOKKOS_INLINE_FUNCTION void update(double d_t, MC::ParticleDataHolder &p, const LocalConcentrationView &concentration, MC::KPRNG _rng)
    {

      const double s = Kokkos::max(0., concentration(0));
      const double phi_s_pts = phi_pts(a_pts, s);
      const double phi_s_in = phi_s_pts + phi_permease(n_permease, a_permease, s);
      const double gamma_PTS_S = phi_s_pts / phi_pts_max;
      const double nu_p = phi_s_in / factor / Yields::x_s;
      const double nu_eff = Kokkos::min(nu_p, nu);
      const double phi_o2_in = phi_o2_max * mass();           // TODO use min(max,phi_phy)
      const double phi_a_in = phi_a_max * MolarMass::acetate; // TODO use min(max,phi_phy)

      metabolism(phi_s_in, phi_o2_in, phi_a_in);
      const double lengthening = q[_INDICES(NU)];
      this->lenght += d_t * lengthening;
      this->nu += d_t * (1.0 / tau_metabolism) * (nu_p - this->nu);
      this->a_pts += d_t * 1.0 / Tau::pts * (MONOD_RATIO(1., s, k_pts) - this->a_pts);
      this->a_permease += d_t * (((1.0 / Tau::Au) * gamma_PTS_S + (1.0 / Tau::Ad) * (1.0 - gamma_PTS_S)) * (1.0 - gamma_PTS_S - a_permease));

      this->n_permease += d_t * (MONOD_RATIO(1. / Tau::new_permease, k_pts, s) + MONOD_RATIO(1. / Tau::tau_rm_perm, s, k_pts)) *
                          (MONOD_RATIO(NPermease_max, k_pts, s) - n_permease);

      std::copy(q.begin(), q.end() - 1, contribs.begin());

      Models::update_division_status(p.status, d_t, GammaDivision::threshold_linear(lenght, l_0, l_1), _rng);
    }

    KOKKOS_INLINE_FUNCTION UptakeAcetate division(MC::ParticleDataHolder &p, MC::KPRNG _rng)
    {
      const double original_lenght = this->lenght;
      const double original_n_permease = this->n_permease;
      this->lenght = original_lenght / 2.;
      this->n_permease = original_n_permease / 2;

      auto child = *this;
      const bool mask = a_pts > a_permease;
      const double max_perm = mask ? a_pts : 1;
      const double max_pts = mask ? 1 : a_permease;
      {
        auto generator = _rng.random_pool.get_state();
        child.a_pts = Kokkos::min(max_pts, Kokkos::max(generator.normal(this->a_pts, this->a_pts / 2.), 0.));
        child.a_permease = Kokkos::min(max_perm, Kokkos::max(generator.normal(this->a_permease, this->a_permease / 2.), 0.));
        _rng.random_pool.free_state(generator);
      }

      return child;
    }

    KOKKOS_INLINE_FUNCTION void contribution(MC::ParticleDataHolder &p, ContributionView contribution)
    {
      auto access_contribs = contribution.access();
      //   access_contribs(0, p.current_container) += -this->contrib_s * p.weight;

      for (int i = 0; i < contribs.size(); ++i)
      {
        access_contribs(i, p.current_container) += contribs[i] * p.weight;
      }
    }

    [[nodiscard]] KOKKOS_INLINE_FUNCTION double mass() const noexcept
    {
      return factor * lenght;
    }

    inline model_properties_detail_t get_properties()
    {

      return {{"lenght", lenght}, {"a_pts", a_pts}, {"mu_eff", nu}, {"a_permease", a_permease}, {"n_permease", n_permease}, {"mass", mass()}};
    }
  };
} // namespace Models

#endif