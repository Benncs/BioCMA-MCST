#ifndef __MODELS_SIMPLE_ECOLI_HPP__
#define __MODELS_SIMPLE_ECOLI_HPP__

#include <mc/particles/particle_model.hpp>
#include <models/uptake.hpp>
#include <models/utils.hpp>

/**
  @brief Contant and definition specific to model Meta  
 */
namespace implMeta
{
  constexpr float l_max = 5e-6;
  constexpr float l_c = 3e-6;
  constexpr float tau_1 = 1000.;
  constexpr float tau_2 = 1000.;

  constexpr float a_11 = 1 / 180.;
  constexpr float a_21 = 1 / 32.;

  constexpr float a_12 = 1 / 180.;
  constexpr float a_22 = 1 / 32.;
  constexpr float d_m = 0.6e-6;
  constexpr float minimal_length = 0.9e-6;
  // static_assert(minimal_length>d_m*4, "Diameter too large");
  constexpr float factor = 1000. * 3.14 * d_m * d_m / 4.;
  constexpr float y_sx = 1. / 2.217737e+00; // y1b

  float consteval get_phi_s_max(float nu)
  {
    return (nu * factor) * y_sx; // kgS/s
  }

  constexpr float mt = factor * minimal_length;

  constexpr float mu_max = 0.8 / 3600.;

  constexpr float nu_max2 = mu_max * mt;

  constexpr float nu_max = 8 * 2e-10; // m/s  https://doi.org/10.7554/eLife.67495;

  constexpr float mu_max_ = nu_max * factor / mt * 3600;

  constexpr float phi_s_max = get_phi_s_max(nu_max); // 5*mt/7200;
  constexpr float phi_perm_max = phi_s_max / 40.;    // 5*mt/7200;
  constexpr float phi_o2_max = 10 * phi_s_max / 180 * 3 * 32;
  constexpr float NPermease_max = 200;

  constexpr float y_sxf = y_sx / 3.; // y1
  // constexpr float y_os = 4.432918e-01;
  constexpr float y_sa = 0.8;

  constexpr float y_os = 3; // 3 mol o2 per mol for glucose

  constexpr float k_pts = 1e-3;
  constexpr float k_o = 0.0001; // g/L: Anane et. al 2017 (Biochem. Eng. J)

  constexpr float mmx = 113.1;

  constexpr float MolarMassG = 180;
  constexpr float MolarMassO2 = 32;

  constexpr float alpha_divison = 5.2;
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
      constexpr auto l_c = implMeta::l_c;
      constexpr double minimal_length = implMeta::minimal_length;
      constexpr double l_max = implMeta::l_max;

      nu1 = implMeta::nu_max * implMeta::factor / 4;
      nu2 = 0; // pimpl.nu1/5.;

      contrib = {0., 0.0, 0.};
      auto g = _rng.random_pool.get_state();
      Uptake::distribute_init<float>(*this, g);
      length = static_cast<float>(Kokkos::max(minimal_length, g.normal(minimal_length, 0.7e-6)));
      l_cp = static_cast<float>(
          Kokkos::min(Kokkos::max(minimal_length, g.normal(l_c, l_c / 7.)), l_max));
      _rng.random_pool.free_state(g);
    }

    KOKKOS_FUNCTION void update(double d_t,
                                MC::ParticleDataHolder& p,
                                const LocalConcentrationView& concentrations,
                                MC::KPRNG _rng)
    {
      using namespace implMeta;
      const auto s = static_cast<float>(concentrations(0));

      const auto phi_s = Models::Uptake::uptake(
          static_cast<float>(d_t), *this, s, phi_s_max, phi_perm_max, NPermease_max);

      const auto o = Kokkos::max(static_cast<float>(concentrations(1)), 0.F);
      const float phi_o2 = (phi_o2_max)*o / (o + k_o); // gO2/s
      const float nu_1_star =
          y_sx * MolarMassG * Kokkos::min(phi_s / MolarMassG, phi_o2 / MolarMassO2 / y_os); // gX/s

      const float s_1_star = (1 / y_sx * nu_1_star);
      const float phi_s_residual_1_star = Kokkos::max(phi_s - s_1_star, 0.F);
      KOKKOS_ASSERT(phi_s_residual_1_star >= 0.F);
      const float nu_2_star = y_sxf * phi_s_residual_1_star; // gX/s

      nu_eff_1 = Kokkos::min(nu_1_star, nu1); // gX/s
      const float s_1 = (1 / y_sx * nu_eff_1);
      const float phi_s_residual_1 = Kokkos::max(phi_s - s_1, 0.F);

      nu_eff_2 = Kokkos::min(y_sxf * phi_s_residual_1, nu2); // gX/s

      const float s_growth = s_1 + (1 / y_sxf * nu_eff_2);
      const float s_overflow = phi_s - s_growth;

      KOKKOS_ASSERT(nu_eff_1 >= 0.F);
      KOKKOS_ASSERT(nu_eff_2 >= 0.F);

      p.status = (length > l_cp) ? MC::CellStatus::CYTOKINESIS : MC::CellStatus::IDLE;

      contrib.phi_s = -phi_s;
      contrib.phi_o =
          -1 * ((1. / y_sx / MolarMassG * y_os * MolarMassO2 * nu_eff_1) + 0. * nu_eff_2);
      contrib.phi_a = nu_eff_2 / y_sxf * y_sa + (s_overflow > 0. ? y_sa * (s_overflow) : 0);

      nu1 = nu1 + static_cast<float>(d_t) * ((nu_1_star - nu1) / tau_1);
      nu2 = nu2 + static_cast<float>(d_t) * ((nu_2_star - nu2) / tau_2);
      length = length + static_cast<float>(d_t) * ((nu_eff_1 + nu_eff_2) / factor);
    }

    KOKKOS_FUNCTION Meta2 division(MC::ParticleDataHolder& p, MC::KPRNG _rng)
    {
      (void)p;
      constexpr auto l_c = implMeta::l_c;
      constexpr double minimal_length = implMeta::minimal_length;
      constexpr double l_max = implMeta::l_max;

      const float l = length / 2.F;

      auto child_pimpl = Models::Meta2(*this); // NOLINT

      length = l;
      child_pimpl.length = l;
      auto g = _rng.random_pool.get_state();
      child_pimpl.l_cp = static_cast<float>(
          Kokkos::min(Kokkos::max(minimal_length, g.normal(l_c, l_c / 6.)), l_max));

      Uptake::distribute_division<float>(*this, child_pimpl, g);
      _rng.random_pool.free_state(g);
      return child_pimpl;
    }

    KOKKOS_FUNCTION void contribution(MC::ParticleDataHolder& p,
                                      const ContributionView& contribution)
    {
      auto access_contribs = contribution.access();

      access_contribs(0, p.current_container) += p.weight * contrib.phi_s;
      access_contribs(1, p.current_container) += p.weight * contrib.phi_o;
      access_contribs(2, p.current_container) += p.weight * contrib.phi_a;
    }

    [[nodiscard]] KOKKOS_FUNCTION double mass() const noexcept
    {

      return implMeta::factor * length;
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
      full(1) = length * 1e6;
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