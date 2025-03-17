#ifndef __TWOMETA_MODEL_HPP__
#define __TWOMETA_MODEL_HPP__

#include "common/traits.hpp"
#include "models/utils.hpp"
#include <mc/prng/prng_extension.hpp>
#include <mc/traits.hpp>
#include <models/uptake.hpp>
#include <string_view>

namespace
{
  template <FloatingPointType F> static F consteval get_phi_s_max(F density, F dl)
  {
    // dl and density must be same unit, dl*density -> mass and y is mass yield
    return (dl * density) * 0.5;
  }
} // namespace

namespace Models
{

  struct TwoMetaModel
  {
    using uniform_weight = std::true_type; // Using type alias
    using Self = TwoMetaModel;
    using FloatType = float;

    enum class particle_var : int
    {
      age = INDEX_FROM_ENUM(Uptakeparticle_var::COUNT),
      length,
      nu1,
      nu2,
      l_cp,
      nu_eff_1, // This is not itself a model property but stored to be exported
      nu_eff_2, // This is not itself a model property but stored to be exported
      // TODO FIND BETTER WAY TO STORE/GET CONTRIBUTIONS
      contrib_phi_s,
      contrib_phi_o2,
      contrib_phi_ac,
      COUNT
    };

    static constexpr std::size_t n_var = INDEX_FROM_ENUM(particle_var::COUNT);
    static constexpr std::string_view name = "simple";
    using SelfParticle = MC::ParticlesModel<Self::n_var, Self::FloatType>;

    MODEL_CONSTANT FloatType l_max_m = 5e-6;   // m
    MODEL_CONSTANT FloatType l_c_m = 3e-6;     // m
    MODEL_CONSTANT FloatType d_m = 0.6e-6;     // m
    MODEL_CONSTANT FloatType l_min_m = 0.9e-6; // m
    MODEL_CONSTANT FloatType lin_density = c_linear_density(static_cast<FloatType>(1000), d_m);
    MODEL_CONSTANT FloatType MolarMassG = Models::MolarMass::GramPerMole::glucose<float>;
    MODEL_CONSTANT FloatType MolarMassO2 = Models::MolarMass::GramPerMole::dioxygen<float>; // g/mol

    MODEL_CONSTANT FloatType y_sx_1 = 1. / 2.217737e+00; // Mode 1 S to X yield (mass)
    MODEL_CONSTANT FloatType y_sx_2 = y_sx_1 / 3.;       // Mode 2 S to X yield (mass)
    MODEL_CONSTANT FloatType y_sa = 0.8;                 // S to A yield (mass)
    MODEL_CONSTANT FloatType y_os_molar = 3;             // 3 mol o2 per mol for glucose
    MODEL_CONSTANT FloatType k_o = 0.0001; // g/L: Anane et. al 2017 (Biochem. Eng. J) (g/g)
    MODEL_CONSTANT FloatType dl_max_ms = 8 * 2e-10; // m/s  https://doi.org/10.7554/eLife.67495;
    MODEL_CONSTANT FloatType tau_1 = 1000.;         // s
    MODEL_CONSTANT FloatType tau_2 = 1000.;         // s

    MODEL_CONSTANT FloatType phi_s_max = (dl_max_ms * lin_density) * y_sx_1;
    MODEL_CONSTANT FloatType phi_perm_max = phi_s_max / 40.; // kgS/
    MODEL_CONSTANT FloatType phi_o2_max =
        10 * phi_s_max / MolarMassG * y_os_molar * MolarMassO2; // kgS/s

    MODEL_CONSTANT FloatType nu_max_kg_s = dl_max_ms * lin_density;

    KOKKOS_INLINE_FUNCTION static void
    init(const MC::KPRNG::pool_type& random_pool, std::size_t idx, const SelfParticle& arr);

    KOKKOS_INLINE_FUNCTION static MC::Status update(const MC::KPRNG::pool_type& random_pool,
                                                    FloatType d_t,
                                                    std::size_t idx,
                                                    const SelfParticle& arr,
                                                    const MC::LocalConcentration& c);

    KOKKOS_INLINE_FUNCTION static void division(const MC::KPRNG::pool_type& random_pool,
                                                std::size_t idx,
                                                std::size_t idx2,
                                                const SelfParticle& arr,
                                                const SelfParticle& buffer_arr);

    KOKKOS_INLINE_FUNCTION static void contribution(std::size_t idx,
                                                    std::size_t position,
                                                    double weight,
                                                    const SelfParticle& arr,
                                                    const MC::ContributionView& contributions);

    KOKKOS_INLINE_FUNCTION static double mass(std::size_t idx, const SelfParticle& arr)
    {
      return GET_PROPERTY(Self::particle_var::length) * lin_density;
    }

    inline constexpr static std::array<std::string_view, n_var> names()
    {
      constexpr std::size_t ln_var = n_var - Uptake<Self>::n_var;
      constexpr auto _names = concat_arrays<Uptake<Self>::n_var, ln_var>(Uptake<Self>::names(),
                                                                         {
                                                                             "age",
                                                                             "length",
                                                                             "nu1",
                                                                             "nu2",
                                                                             "l_cp",
                                                                             "nu_eff_1",
                                                                             "nu_eff_2",
                                                                             "contrib_phi_s",
                                                                             "contrib_phi_o2",
                                                                             "contrib_phi_ac",
                                                                         });

      return _names;
    }
  };

  CHECK_MODEL(TwoMetaModel)

  KOKKOS_INLINE_FUNCTION void
  TwoMetaModel::init([[maybe_unused]] const MC::KPRNG::pool_type& random_pool,
                     std::size_t idx,
                     const SelfParticle& arr)
  {

    constexpr auto length_dist = MC::Distributions::TruncatedNormal<FloatType>(
        l_min_m, l_min_m / 5., 0.5 * l_min_m, l_max_m);

    constexpr auto length_c_dist =
        MC::Distributions::TruncatedNormal<FloatType>(l_c_m, l_c_m / 7., l_min_m, l_max_m);

    auto gen = random_pool.get_state();
    GET_PROPERTY(Self::particle_var::length) = length_dist.draw(gen);
    GET_PROPERTY(Self::particle_var::l_cp) = length_c_dist.draw(gen);
    random_pool.free_state(gen);
    GET_PROPERTY(particle_var::contrib_phi_s) = 0;
    Uptake<Self>::init(random_pool, idx, arr);
  }

  KOKKOS_INLINE_FUNCTION MC::Status
  TwoMetaModel::update(const MC::KPRNG::pool_type& random_pool,
                       FloatType d_t,
                       std::size_t idx,
                       const SelfParticle& arr,
                       const MC::LocalConcentration& concentrations)
  {
    (void)random_pool;
    const auto phi_s =
        Uptake<Self>::uptake_step(phi_s_max, phi_perm_max, d_t, idx, arr, concentrations);

    GET_PROPERTY(Self::particle_var::age) += d_t;
    GET_PROPERTY(Self::particle_var::contrib_phi_s) = phi_s;

    const auto s = static_cast<float>(concentrations(0));
    const auto o = 9e-3; // Kokkos::max(static_cast<float>(concentrations(1)), 0.F);

    const float phi_o2 = (phi_o2_max)*o / (o + k_o); // gO2/s

    const float nu_1_star = y_sx_1 * MolarMassG *
                            Kokkos::min(phi_s / MolarMassG,
                                        phi_o2 / MolarMassO2 / y_os_molar); // gX/s

    const float s_1_star = (1 / y_sx_1 * nu_1_star);

    const float phi_s_residual_1_star = Kokkos::max(phi_s - s_1_star, 0.F);
    KOKKOS_ASSERT(phi_s_residual_1_star >= 0.F);

    const float nu_2_star = y_sx_2 * phi_s_residual_1_star; // gX/s
    GET_PROPERTY(Self::particle_var::nu_eff_1) =
        Kokkos::min(nu_1_star, GET_PROPERTY(Self::particle_var::nu1)); // gX/s
    const float s_1 = (1 / y_sx_1 * GET_PROPERTY(Self::particle_var::nu_eff_1));
    const float phi_s_residual_1 = Kokkos::max(phi_s - s_1, 0.F);
    GET_PROPERTY(Self::particle_var::nu_eff_2) =
        Kokkos::min(y_sx_2 * phi_s_residual_1, GET_PROPERTY(Self::particle_var::nu2)); // gX/s

    const float s_growth = s_1 + (1 / y_sx_2 * GET_PROPERTY(Self::particle_var::nu_eff_2));

    const float s_overflow = phi_s - s_growth;

    KOKKOS_ASSERT(GET_PROPERTY(Self::particle_var::nu_eff_1) >= 0.F);
    KOKKOS_ASSERT(GET_PROPERTY(Self::particle_var::nu_eff_2) >= 0.F);

    GET_PROPERTY(Self::particle_var::contrib_phi_s) = -phi_s;
    GET_PROPERTY(Self::particle_var::contrib_phi_o2) =
        -1 * ((1. / y_sx_1 / MolarMassG * y_os_molar * MolarMassO2 *
               GET_PROPERTY(Self::particle_var::nu_eff_1)) +
              0. * GET_PROPERTY(Self::particle_var::nu_eff_2));
    GET_PROPERTY(Self::particle_var::contrib_phi_ac) =
        GET_PROPERTY(Self::particle_var::nu_eff_2) / y_sx_2 * y_sa +
        (s_overflow > 0. ? y_sa * (s_overflow) : 0);

    GET_PROPERTY(Self::particle_var::nu1) =
        GET_PROPERTY(Self::particle_var::nu1) +
        static_cast<float>(d_t) * ((nu_1_star - GET_PROPERTY(Self::particle_var::nu1)) / tau_1);

    GET_PROPERTY(Self::particle_var::nu2) =
        GET_PROPERTY(Self::particle_var::nu2) +
        static_cast<float>(d_t) * ((nu_2_star - GET_PROPERTY(Self::particle_var::nu2)) / tau_2);

    GET_PROPERTY(Self::particle_var::length) =
        GET_PROPERTY(Self::particle_var::length) +
        static_cast<float>(d_t) * ((GET_PROPERTY(Self::particle_var::nu_eff_1) +
                                    GET_PROPERTY(Self::particle_var::nu_eff_2)) /
                                   (lin_density));

    return (GET_PROPERTY(Self::particle_var::length) > GET_PROPERTY(Self::particle_var::l_cp))
               ? MC::Status::Division
               : MC::Status::Idle;
  }

  KOKKOS_INLINE_FUNCTION void TwoMetaModel::division(const MC::KPRNG::pool_type& random_pool,
                                                     std::size_t idx,
                                                     std::size_t idx2,
                                                     const SelfParticle& arr,
                                                     const SelfParticle& child_buffer_arr)
  {

    const FloatType new_current_length =
        GET_PROPERTY(particle_var::length) / static_cast<FloatType>(2.);

    GET_PROPERTY(Self::particle_var::length) = new_current_length;
    GET_PROPERTY(Self::particle_var::age) = 0;

    GET_PROPERTY_FROM(idx2, child_buffer_arr, Self::particle_var::length) = new_current_length;
    GET_PROPERTY_FROM(idx2, child_buffer_arr, Self::particle_var::age) = 0;

    constexpr auto length_c_dist =
        MC::Distributions::TruncatedNormal<FloatType>(l_c_m, l_c_m / 7., l_min_m, l_max_m);

    auto gen = random_pool.get_state();
    GET_PROPERTY(Self::particle_var::l_cp) = length_c_dist.draw(gen);
    GET_PROPERTY_FROM(idx2, child_buffer_arr, Self::particle_var::l_cp) = length_c_dist.draw(gen);
    random_pool.free_state(gen);

    Uptake<Self>::division(random_pool, idx, idx2, arr, child_buffer_arr);
  }

  KOKKOS_INLINE_FUNCTION void TwoMetaModel::contribution([[maybe_unused]] std::size_t idx,
                                                         std::size_t position,
                                                         double weight,
                                                         [[maybe_unused]] const SelfParticle& arr,
                                                         const MC::ContributionView& contributions)
  {
    auto access = contributions.access();
    access(position, 0) += weight * GET_PROPERTY(Self::particle_var::contrib_phi_s); // NOLINT
  }

} // namespace Models

#endif