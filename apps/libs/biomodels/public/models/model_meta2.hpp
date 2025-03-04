#ifndef __MODELS_SIMPLE_ECOLI_HPP__
#define __MODELS_SIMPLE_ECOLI_HPP__

#include <mc/particles/particle_model.hpp>
#include <mc/prng/prng_extension.hpp>
#include <models/uptake.hpp>
#include <models/utils.hpp>

/**
  @brief Contant and definition specific to model Meta
 */
namespace implMeta
{
  constexpr float m_to_micron = 1e6;
  constexpr float micron_to_m = 1e-6;
  constexpr float g_to_pico = 1e12;
  constexpr float kg_to_pico = g_to_pico * 1e3; // pg->g 1e12 g->kg 1e3
  constexpr float pico_to_kg = 1. / kg_to_pico; // pg->g 1e12 g->kg 1e3

  constexpr float l_max_m = 5e-6;   // m
  constexpr float l_c_m = 3e-6;     // m
  constexpr float d_m = 0.6e-6;     // m
  constexpr float l_min_m = 0.9e-6; // m
  constexpr float MolarMassG = Models::MolarMass::GramPerMole::glucose<float>;
  constexpr float MolarMassO2 = Models::MolarMass::GramPerMole::dioxygen<float>;      // g/mol
  constexpr float dl_max_ms = 8 * 2e-10; // m/s  https://doi.org/10.7554/eLife.67495;
  constexpr float tau_1 = 1000.;         // s
  constexpr float tau_2 = 1000.;         // s
  constexpr float cell_density = 1000.;  // kg/m3
  constexpr float linear_density_kg_m =
      Models::c_linear_density(cell_density, d_m); // kg/m We stay in meter
  constexpr float NPermease_max = 200;             // Number of permeases
  // Yields
  constexpr float y_sx_1 = 1. / 2.217737e+00; // Mode 1 S to X yield (mass)
  constexpr float y_sx_2 = y_sx_1 / 3.;       // Mode 2 S to X yield (mass)
  constexpr float y_sa = 0.8;                 // S to A yield (mass)
  constexpr float y_os_molar = 3;             // 3 mol o2 per mol for glucose
  constexpr float k_o = 0.0001;               // g/L: Anane et. al 2017 (Biochem. Eng. J) (g/g)

  // Uptake
  float consteval get_phi_s_max(float density, float dl)
  {
    // dl and density must be same unit, dl*density -> mass and y is mass yield
    return (dl * density) * y_sx_1;
  }

  constexpr float phi_s_max = get_phi_s_max(linear_density_kg_m, dl_max_ms); // kgS/s
  constexpr float phi_perm_max = phi_s_max / 40.;                            // kgS/s
  constexpr float phi_o2_max = 10 * phi_s_max / MolarMassG *
                               y_os_molar * MolarMassO2;                   // kgS/s
  constexpr float nu_max_kg_s = dl_max_ms * implMeta::linear_density_kg_m; // kg/s

  // Conversions
  constexpr float l_max_um = l_max_m * m_to_micron;
  constexpr float l_c_um = l_c_m * m_to_micron;
  constexpr float d_m_um = d_m * m_to_micron;
  constexpr float l_min_um = l_min_m * m_to_micron;
  constexpr float linear_density_kg_um = linear_density_kg_m * micron_to_m; // kg/m
  constexpr float dl_max_um_s = dl_max_ms * m_to_micron;                    // um/s
  constexpr float nu_max_pg_s = nu_max_kg_s * kg_to_pico;
  constexpr float phi_s_max_pg = phi_s_max * kg_to_pico;       // pgS/s
  constexpr float phi_perm_max_pg = phi_perm_max * kg_to_pico; // pgS/s
  constexpr float phi_o2_max_pg = phi_o2_max * kg_to_pico;     // pgS/s

  constexpr auto length_dist =
      MC::Distributions::TruncatedNormal<float>(l_min_um, l_min_um / 5., 0.5 * l_min_um, l_max_um);

  constexpr auto length_c_dist =
      MC::Distributions::TruncatedNormal<float>(l_c_um, l_c_um / 7., l_min_um, l_max_um);

  constexpr auto mu_nu_dist = nu_max_pg_s * 0.1;
  constexpr float nu_dist_factor = 1e4; // Scale factor
  constexpr auto nu_1_initial_dist = MC::Distributions::TruncatedNormal<float>(
      mu_nu_dist, mu_nu_dist / 7., 0., static_cast<double>(nu_max_pg_s));

} // namespace implMeta

namespace Models
{

  /**
  @brief Two growth rate metabolic model
   */
  struct Meta2
  {

    struct contribs
    {
      float phi_s;
      float phi_o;
      float phi_a;
    };

    DECLARE_UPTAKE_PROPERTIES(float)
    float length;
    float nu1;
    float nu2;
    float l_cp;
    float mu;
    float nu_eff_1;
    float nu_eff_2;
    contribs contrib;

    KOKKOS_FUNCTION void init(MC::ParticleDataHolder& p, MC::KPRNG _rng)
    {
      (void)p;

      constexpr auto length_dist = implMeta::length_dist;
      constexpr auto length_c_dist = implMeta::length_c_dist;
      constexpr auto nu_dist = implMeta::nu_1_initial_dist;
      nu2 = 0;

      contrib = {0., 0.0, 0.};
      auto g = _rng.random_pool.get_state();
      nu1 = nu_dist.draw(g);
      Uptake::distribute_init<float>(*this, g);
      length = length_dist.draw(g);
      l_cp = length_c_dist.draw(g);
      _rng.random_pool.free_state(g);
    }

    KOKKOS_FUNCTION void update(double d_t,
                                MC::ParticleDataHolder& p,
                                const LocalConcentrationView& concentrations,
                                MC::KPRNG _rng)
    {
      (void)_rng;
      using namespace implMeta;
      const auto s = static_cast<float>(concentrations(0));

      const auto phi_s = Models::Uptake::uptake(
          static_cast<float>(d_t), *this, s, phi_s_max_pg, phi_perm_max_pg, NPermease_max);

      const auto o = Kokkos::max(static_cast<float>(concentrations(1)), 0.F);

      const float phi_o2 = (phi_o2_max_pg)*o / (o + k_o); // gO2/s

      const float nu_1_star = y_sx_1 * MolarMassG*
                              Kokkos::min(phi_s / MolarMassG,
                                          phi_o2 / MolarMassO2 / y_os_molar); // gX/s

      const float s_1_star = (1 / y_sx_1 * nu_1_star);

      const float phi_s_residual_1_star = Kokkos::max(phi_s - s_1_star, 0.F);
      KOKKOS_ASSERT(phi_s_residual_1_star >= 0.F);

      const float nu_2_star = y_sx_2 * phi_s_residual_1_star; // gX/s
      nu_eff_1 = Kokkos::min(nu_1_star, nu1);                 // gX/s
      const float s_1 = (1 / y_sx_1 * nu_eff_1);
      const float phi_s_residual_1 = Kokkos::max(phi_s - s_1, 0.F);
      nu_eff_2 = Kokkos::min(y_sx_2 * phi_s_residual_1, nu2); // gX/s

      const float s_growth = s_1 + (1 / y_sx_2 * nu_eff_2);

      const float s_overflow = phi_s - s_growth;

      KOKKOS_ASSERT(nu_eff_1 >= 0.F);
      KOKKOS_ASSERT(nu_eff_2 >= 0.F);

      p.status = (length > l_cp) ? MC::CellStatus::CYTOKINESIS : MC::CellStatus::IDLE;

      contrib.phi_s = -phi_s;
      contrib.phi_o = -1 * ((1. / y_sx_1 / MolarMassG *
                             y_os_molar * MolarMassO2 * nu_eff_1) +
                            0. * nu_eff_2);
      contrib.phi_a = nu_eff_2 / y_sx_2 * y_sa + (s_overflow > 0. ? y_sa * (s_overflow) : 0);

      nu1 = nu1 + static_cast<float>(d_t) * ((nu_1_star - nu1) / tau_1);
      nu2 = nu2 + static_cast<float>(d_t) * ((nu_2_star - nu2) / tau_2);
      length = length + static_cast<float>(d_t) *
                            ((nu_eff_1 + nu_eff_2) / (linear_density_kg_um * kg_to_pico));
    }

    KOKKOS_FUNCTION Meta2 division(MC::ParticleDataHolder& p, MC::KPRNG _rng)
    {
      (void)p;

      const float l = 0.5F * length;

      constexpr auto nu_dist_factor = implMeta::nu_dist_factor;

      constexpr auto length_c_dist = implMeta::length_c_dist;
      auto child = Models::Meta2(*this);

      length = l;
      child.length = l;
      auto g = _rng.random_pool.get_state();

      // child.l_cp = length_c_dist.draw(g);

      // child.nu1 = static_cast<float>(MC::Distributions::ScaledTruncatedNormal<double>::draw_from(
      //     g, nu_dist_factor, this->nu1, this->nu1 / 3., 0., 1.));

      const float sigma_nu2 = (almost_equal(nu2, 0., 1e-13)) ? 0.01F : this->nu2 / 3.;

      // child.nu2 =
      // static_cast<float>(
      // MC::Distributions::ScaledTruncatedNormal<double>::draw_from(
      //     g, nu_dist_factor, this->nu2, sigma_nu2, 0., 1.F));

      KOKKOS_ASSERT(child.nu1 >= 0.F);
      KOKKOS_ASSERT(child.nu2 >= 0.F);

      Uptake::distribute_division<float>(*this, child, g);
      _rng.random_pool.free_state(g);
      return child;
    }

    KOKKOS_FUNCTION void contribution(MC::ParticleDataHolder& p,
                                      const ContributionView& contribution)
    {
      auto access_contribs = contribution.access();

      access_contribs(0, p.current_container) += p.weight * contrib.phi_s * implMeta::pico_to_kg;
      access_contribs(1, p.current_container) += p.weight * contrib.phi_o * implMeta::pico_to_kg;
      ;
      access_contribs(2, p.current_container) += p.weight * contrib.phi_a * implMeta::pico_to_kg;
      ;
    }

    [[nodiscard]] KOKKOS_FUNCTION double mass() const noexcept
    {
      // Return the mass in kg
      return implMeta::linear_density_kg_um * length;
    }

    template <class Archive> void serialize(Archive& ar)
    {
      ar(length, nu1, nu2, a_pts, a_permease);
    }

    static std::vector<std::string> names()
    {
      return {"mass", "length", "nu1", "nu2", "nu_eff_1", "nu_eff_2", "a_pts", "a_permease"};
    }

    KOKKOS_INLINE_FUNCTION static std::size_t get_number()
    {
      return 8;
    }

    KOKKOS_INLINE_FUNCTION void fill_properties(SubViewtype full) const
    {
      full(0) = mass();
      full(1) = length;
      full(2) = nu1;
      full(3) = nu2;
      full(4) = nu_eff_1;
      full(5) = nu_eff_2;
      full(6) = a_pts;
      full(7) = a_permease;
    }
  };
} // namespace Models

#endif