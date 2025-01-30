#ifndef __BIO_MODEL_UPTAKE_ACETATE_HPP__
#define __BIO_MODEL_UPTAKE_ACETATE_HPP__

#include "mc/prng/prng.hpp"
#include <Kokkos_Core.hpp>
#include <common/kokkos_vector.hpp>
#include <mc/particles/particle_model.hpp>

#include <models/utils.hpp>
#include <numbers>

namespace implEcoli
{

#define MONOD_RATIO(__c1__, __x1__, __k1__) ((__c1__) * (__x1__) / ((__x1__) + (__k1__)))

#define _INDICES(__name__) static_cast<std::size_t>(Indices::__name__)

  struct MolarMass
  {
    static constexpr double glucose = 180e-3;
    static constexpr double dioxygen = 32e-3;
    static constexpr double acetate = 59e-3;
    static constexpr double co2 = 44e-3;
    static constexpr double X = 113.1e-3;
  };
  struct Yields
  {

    static constexpr double Y_EG = 12.05; // mol_E/mol_G
    static constexpr double Y_XG = 1.15;  // mol_X/mol_G
    static constexpr double m = 250;      // µmol_G/g_X.h
    static constexpr double Y_OG = 6.;    // mol_O/mol_G
    static constexpr double Yo_EG = 20.;  // mol_E/mol_G
    static constexpr double Yf_EG = 3.;   // mol_E/mol_G
    static constexpr double Y_AG = 1.5;   // mol_A/mol_G (modified for CO2 production)
    // static constexpr double Y_CO2G_over = 2.; // mol_CO2/mol_G (overflow metabolism) not used

    static constexpr double Y_EA = 12.05 / 3.; // mol_E/mol_A
    static constexpr double Y_XA = 0.4;        // mol_X/mol_A
    static constexpr double Y_OA = 2.;         // mol_O/mol_A
    static constexpr double Yo_EA = 4.67;      // mol_E/mol_A
    static constexpr double x_s = 1.15;        // 2.;

    static constexpr double nu_s_x = MolarMass::glucose / MolarMass::X *
                                     (Yields::Y_EG + Yields::Yo_EG) /
                                     (Yields::Y_XG * Yields::Yo_EG);
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

  constexpr double l_1 = 4e-6;
  constexpr double l_0 = 3e-6;
  constexpr double d_m = l_0 / 5.;
  constexpr double minimal_length = d_m * 1.4;
  constexpr double factor = 1000. * 3.14 * d_m * d_m / 4.;
  static_assert(l_0 < l_1, "Length");

  consteval double f_num_max()
  {
    return l_0 * std::numbers::ln2 / 3600.; // NOLINT
  }

  consteval double f_phi_pts_max(double nu_m)
  {
    return nu_m * Yields::x_s * factor;
  }

  constexpr double nu_max = f_num_max(); // m/s

  constexpr double phi_pts_max = f_phi_pts_max(nu_max);
  constexpr double phi_o2_max = 15.6 / 3600;     // mol0/gx/s
  constexpr double phi_a_max = phi_pts_max / 3.; // mola/kgx/s
  constexpr double k_pts = 1e-3;
  constexpr double kppermease = 1e-2;
  constexpr double phi_permease_specific = phi_pts_max / 40.;

  constexpr double psi_o_meta = 20e-3 / 3600 * 32e-3; // 20mmmol O2 /h;
  constexpr double tau_metabolism = 3600;

  constexpr double NPermease_max = 200;
  constexpr double NPermease_init = 1.;

  KOKKOS_INLINE_FUNCTION double phi_pts(double a_pts, double S)
  {
    return a_pts * MONOD_RATIO(phi_pts_max, S, k_pts);
  }

  KOKKOS_INLINE_FUNCTION double phi_permease(double n_permease, double a_permease, double S)
  {
    return a_permease * MONOD_RATIO(n_permease * phi_permease_specific, S, kppermease);
  }
  // static_assert(NPermease_max * phi_permease_specific == 5. * phi_pts_max, "pts/perm");

} // namespace implEcoli

namespace Models
{

  struct Ecoli
  {

    using PhiUptakes = struct Contribs
    {
      double glucose;
      double acetate;
      double oxygen;
    };

    struct Rates
    {
      double glucose;
      double acetate;
      double oxygen;
      double carbon_dioxide;
      double nu;
    };

    static constexpr std::size_t n_meta = 5;

    enum class Indices : std::size_t
    {
      NU = 4,
      GLUCOSE = 0,
      O2 = 1,
      Ac = 2,
      CO2 = 3
    };

    double length;
    double nu;
    double a_pts;
    double a_permease;
    double n_permease;
    Rates rates;
    PhiUptakes phi_uptakes;

    KOKKOS_FUNCTION void init(MC::ParticleDataHolder& p, MC::KPRNG _rng)
    {
      using namespace implEcoli;
      constexpr double local_l = minimal_length;
      auto generator = _rng.random_pool.get_state();

      this->length =
          Kokkos::max(local_l, Kokkos::max(generator.normal(l_0 / 1.5, l_0 / 1.5 / 5.), 0.));
      this->a_permease = Kokkos::max(generator.normal(1e-3, 1e-4), 0.);
      this->a_pts = Kokkos::min(1., Kokkos::max(generator.normal(0.8, 0.1), 0.));
      this->n_permease =
          Kokkos::max(generator.normal(NPermease_init / 2., NPermease_init / 5.), 0.);

      this->nu =
          nu_max; // Kokkos::max(generator.normal(nu_max / 2., nu_max / 2. / 5.), 0.) / length;
      _rng.random_pool.free_state(generator);
    }

    KOKKOS_FUNCTION void metabolism()
    {
      using namespace ::implEcoli;
      // yield n for mode i, yields in mol/mol
      // E: ATP/G:glucose/O:Oxygen/A:Acetate
      constexpr double y1_b = (Yields::Y_EG + Yields::Yo_EG) / Yields::Yo_EG; // None
      constexpr double y2_b =
          (Yields::Y_EG + Yields::Yo_EG) /
          (Yields::Y_EG * Yields::Y_OG); // G/O //To modify divide by YO_EG not Y_Eg
      constexpr double y1_c = (Yields::Y_EG + Yields::Yf_EG) / Yields::Yf_EG; // None
      constexpr double y1_e = (Yields::Y_EA + Yields::Yo_EA) / Yields::Yo_EA; // None
      constexpr double y2_e =
          (Yields::Y_EA + Yields::Yo_EA) / (Yields::Yo_EA * Yields::Y_OA); // A/O
      constexpr double yield_mol_S_O =
          1.; // Number of mol of S from mole of O for mode D conversion

      // Yield for specie for specifc mode: id oxygen mode b => y_o2_b
      constexpr double y_o2_b = (Yields::Y_EG * Yields::Y_OG) / (Yields::Y_EG + Yields::Yo_EG);
      constexpr double y_o2_e = (Yields::Y_EA * Yields::Y_OA) / (Yields::Y_EA + Yields::Yo_EA);
      constexpr double y_ac_c = (Yields::Y_AG * Yields::Y_EG) / (Yields::Y_EG + Yields::Yf_EG);
      constexpr double y_x_b =
          Yields::Y_XG * Yields::Yo_EG / (Yields::Y_EG + Yields::Yo_EG); // Y_xg /y1_b
      constexpr double y_x_c = Yields::Y_XG * Yields::Yf_EG / (Yields::Y_EG + Yields::Yf_EG);
      constexpr double y_x_e = Yields::Y_XA * Yields::Yo_EA / (Yields::Y_EA + Yields::Yo_EA);

      constexpr double y_phi_o2_residual =
          (Yields::Y_EG * Yields::Y_OG) / (Yields::Y_EG + Yields::Yo_EG); // y2_bthe same

      // Incoming for uptake mecanism
      const auto [phi_s_in, phi_a_in, phi_o2_in] = phi_uptakes; // kg/s
      const double current_mass = mass();

      // Incomin in molee y2_b = (Yields::Y_EG + Yields::Yo_EG) / (Yields::Y_EG * Yields::Y_OG);  //
      // G/O
      const double phi_s_in_mol = phi_s_in / (current_mass * MolarMass::glucose); // mol/s
      const double phi_a_in_mol = phi_a_in / (current_mass * MolarMass::acetate);
      const double phi_o2_in_mol = phi_o2_in;

      // Maximum growth on oxidative (best yield)
      const double r1_max = nu / MolarMass::X / Yields::x_s / MolarMass::glucose; // molS/s
      // Oxydative growth

      const double phi_mode_b =
          min_var(y1_b * r1_max, phi_s_in_mol, y2_b * phi_o2_in_mol); // molS/s

      double phi_s_residual = phi_s_in_mol - phi_mode_b; // molS/s

      const double phi_o2_residual = phi_o2_in_mol - (y_phi_o2_residual * phi_mode_b); // molO2/s

      const double nu_residual =
          ((Yields::Y_EG + Yields::Yo_EG) / Yields::Yo_EG * r1_max) - phi_mode_b; // molS/s

      const double phi_mode_c = min_var(y1_c * nu_residual, phi_s_residual); // molS/s

      phi_s_residual = phi_s_residual - phi_mode_c; // molS/s

      const double phi_mode_d = min_var(phi_s_residual, yield_mol_S_O * phi_o2_residual); // molS/s

      const double phi_mode_e = min_var(y1_e * nu_residual, phi_a_in_mol, y2_e * phi_o2_residual);

      rates.glucose = -phi_s_in;
      rates.oxygen = -current_mass * MolarMass::dioxygen *
                     (y_o2_b * phi_mode_b + y_o2_e * phi_mode_e + phi_mode_d);

      rates.acetate = current_mass * MolarMass::acetate *
                      (y_ac_c * phi_mode_c + Yields::Y_AG * phi_mode_d - phi_mode_e);

      rates.carbon_dioxide = -current_mass * MolarMass::co2 * (rates.oxygen / MolarMass::dioxygen);

      rates.nu = MolarMass::X * (y_x_b * phi_mode_b + y_x_c * phi_mode_c + y_x_e * phi_mode_e);
    }

    KOKKOS_FUNCTION double uptake(const LocalConcentrationView& concentration)
    {
      using namespace implEcoli;
      const double s = Kokkos::max(0., concentration(_INDICES(GLUCOSE)));
      const double phi_s_pts = phi_pts(a_pts, s);
      const double gamma_PTS_S = phi_s_pts / phi_pts_max;
      const double micro_mixing = 0.3;

      phi_uptakes.glucose = phi_s_pts + phi_permease(n_permease, a_permease, s);
      phi_uptakes.acetate = phi_a_max * concentration(_INDICES(Ac)) / 1;
      // 1 / micro_mixing * concentration(_INDICES(Ac)) / 1000.;

      phi_uptakes.oxygen = phi_o2_max; //* concentration(_INDICES(O2)) / 9e-3;
                                       // 1 / micro_mixing * concentration(_INDICES(O2)) / 1000.;
      return gamma_PTS_S;
    }

    KOKKOS_FUNCTION void update(double d_t,
                                MC::ParticleDataHolder& p,
                                const LocalConcentrationView& concentration,
                                Kokkos::Random_XorShift64_Pool<> _rng)
    {
      using namespace implEcoli;
      const double s = Kokkos::max(0., concentration(0));
      const double gamma_PTS_S = uptake(concentration);
      const double nu_p = phi_uptakes.glucose / mass() / (Yields::nu_s_x);

      metabolism();
      this->length += d_t * length * rates.nu;
      this->nu += d_t * (1.0 / tau_metabolism) * (nu_p - this->nu);
      this->a_pts += d_t * 1.0 / Tau::pts * (MONOD_RATIO(1., s, k_pts) - this->a_pts);
      this->a_permease +=
          d_t * (((1.0 / Tau::Au) * gamma_PTS_S + (1.0 / Tau::Ad) * (1.0 - gamma_PTS_S)) *
                 (1.0 - gamma_PTS_S - a_permease));

      this->n_permease += d_t *
                          (MONOD_RATIO(1. / Tau::new_permease, k_pts, s) +
                           MONOD_RATIO(1. / Tau::tau_rm_perm, s, k_pts)) *
                          (MONOD_RATIO(NPermease_max, k_pts, s) - n_permease);

      Models::update_division_status(
          p.status, d_t, GammaDivision::threshold_linear(length, l_0, l_1), _rng);
    }

    KOKKOS_FUNCTION Ecoli division(MC::ParticleDataHolder& p, MC::KPRNG _rng)
    {
      using namespace implEcoli;
      const double original_length = this->length;
      const double original_n_permease = this->n_permease;
      this->length = original_length / 2.;
      this->n_permease = original_n_permease / 2;

      auto child = *this;
      const bool mask = a_pts > a_permease;
      const double max_perm = mask ? a_pts : 1;
      const double max_pts = mask ? 1 : a_permease;

      {
        // auto generator = _rng.random_pool.get_state();
        // child.a_pts =
        //     Kokkos::min(max_pts, Kokkos::max(generator.normal(this->a_pts, this->a_pts / 2.),
        //     0.));

        // child.a_permease = Kokkos::min(
        //     max_perm, Kokkos::max(generator.normal(this->a_permease, this->a_permease / 2.),
        //     0.));

        auto generator = _rng.random_pool.get_state();
        child.a_pts =
            Kokkos::min(1., Kokkos::max(generator.normal(this->a_pts, this->a_pts / 2.), 0.));

        child.a_permease = Kokkos::min(
            1., Kokkos::max(generator.normal(this->a_permease, this->a_permease / 2.), 0.));

        child.nu = Kokkos::max(generator.normal(this->nu, this->nu / 3.), 0.);
        _rng.random_pool.free_state(generator);
      }

      return child;
    }

    KOKKOS_FUNCTION void contribution(MC::ParticleDataHolder& p,
                                      const ContributionView& contribution)
    {
      using namespace implEcoli;
      auto access_contribs = contribution.access();
      access_contribs(_INDICES(GLUCOSE), p.current_container) += rates.glucose * p.weight;
      // access_contribs(_INDICES(Ac), p.current_container) += rates.acetate * p.weight;
      // access_contribs(_INDICES(O2), p.current_container) += rates.oxygen * p.weight;
      // access_contribs(_INDICES(CO2), p.current_container) += rates.carbon_dioxide * p.weight;
    }

    [[nodiscard]] KOKKOS_FUNCTION double mass() const noexcept
    {
      using namespace implEcoli;
      return implEcoli::factor * length;
    }

    template <class Archive> void serialize(Archive& ar)
    {
      ar(length, nu, a_pts, a_permease);
    }

    static std::vector<std::string> names()
    {
      return {"mass", "length", "nu_eff", "nu","a_permease","a_pts"};
    }

    KOKKOS_INLINE_FUNCTION static consteval std::size_t get_number()
    {
      return 6;
    }

    KOKKOS_INLINE_FUNCTION void fill_properties(SubViewtype full) const
    {
      full(0) = mass();
      full(1) = this->length;
      full(2) = this->rates.nu;
      full(3) = this->nu;
      full(4)=this->a_permease;
      full(5)=this->a_pts;
    }
  };

} // namespace Models

#endif