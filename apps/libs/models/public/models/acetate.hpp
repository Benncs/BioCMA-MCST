#ifndef __ACETATE_MODEL_HPP__
#define __ACETATE_MODEL_HPP__

#include "Kokkos_Assert.hpp"
#include "Kokkos_Core_fwd.hpp"
#include "Kokkos_Macros.hpp"
#include "common/common.hpp"
#include "common/traits.hpp"
#include "mc/macros.hpp"
#include "models/utils.hpp"
#include <array>
#include <mc/prng/prng_extension.hpp>
#include <mc/traits.hpp>
#include <optional>
#include <string_view>

namespace Models
{

  /**
   * @brief Compile time elemental closure of the reactions of Models::Acetate.
   *
   * Only the growth yields are free parameters; the oxygen and acetate
   * coefficients are the *only* ones that balance carbon and electrons, so
   * they are derived here instead of being tabulated. Everything is per
   * C-mole of the carbon substrate, with NH3 as nitrogen source; water closes
   * H and O automatically once C, N and the electrons balance.
   *
   * Nothing in this namespace is evaluated at run time.
   */
  namespace Stoichiometry
  {
    inline constexpr double m_o2 = 32.0;           ///< [g/mol]
    inline constexpr double m_cmol_glucose = 30.0; ///< C6H12O6, 180/6 [g/Cmol]
    inline constexpr double m_cmol_acetate = 30.0; ///< C2H4O2, 60/2   [g/Cmol]
    inline constexpr double m_cmol_biomass = 24.6; ///< CH1.8O0.5N0.2  [g/Cmol]

    /// Degree of reduction per C-mole
    inline constexpr double gamma_glucose = 4.0; ///< 4 + 2 - 2
    inline constexpr double gamma_acetate = 4.0; ///< (2*4 + 4 - 2*2)/2
    inline constexpr double gamma_biomass = 4.2; ///< 4 + 1.8 - 1.0 - 0.6

    /// Electrons left by 1 g of substrate once y_x g of biomass are built.
    /// They must go somewhere: to O2 if the reaction is oxidative, to acetate
    /// if it is fermentative. This is what forbids an arbitrary acetate yield.
    constexpr double
    e_surplus(double m_cmol_sub, double gamma_sub, double y_x)
    {
      return gamma_sub / m_cmol_sub - gamma_biomass * y_x / m_cmol_biomass;
    }

    /// sub + O2 -> X + CO2 (+H2O) : gO2 per g of substrate.
    /// At y_x = 0 this gives back the complete combustion demand, 192/180 for
    /// glucose: that value is only correct for a reaction that does not grow.
    constexpr double
    y_o2_ox(double m_cmol_sub, double gamma_sub, double y_x)
    {
      return m_o2 * e_surplus(m_cmol_sub, gamma_sub, y_x) / 4.0;
    }

    /// sub -> X + Ac + CO2 (+H2O) : gAc per g of substrate
    constexpr double
    y_ac_fer(double m_cmol_sub, double gamma_sub, double y_x)
    {
      return m_cmol_acetate * e_surplus(m_cmol_sub, gamma_sub, y_x)
             / gamma_acetate;
    }
  } // namespace Stoichiometry

  /**
   * @brief Acetate model with explicit oxygen limitation.
   *
   * Derived from Models::SimpleAcetate, this model adds dissolved oxygen as a
   * third liquid species.
   *
   * @warning The transfer kernel hardcodes the oxygen row of the kla matrix to
   * index 1, so the species ordering (S, O2, Ac) must not be changed.
   *
   * Four lumped reactions, by decreasing priority:
   *
   *  - @c i_go  S + O2 -> X + CO2      respiration
   *  - @c i_ao  Ac + O2 -> X + CO2     acetate re-uptake
   *  - @c i_ov  S -> Ac                overflow, no growth
   *  - @c i_fe  S -> X + Ac + CO2      fermentation
   *
   * @section balance Mass balance
   *
   * Acetate moves according to the two limitations of the cell and to nothing
   * else:
   *
   *  - @b a_p. The glucose that oxygen could have supported but that the
   *    growth capacity cannot handle leaves as acetate through @c i_ov. This
   *    is the classic overflow, active even at oxygen saturation.
   *  - @b O2. The glucose that the growth capacity could have handled but that
   *    oxygen cannot support is fermented through @c i_fe, which also releases
   *    acetate.
   *  - Acetate is taken back up by @c i_ao alone, on the oxygen *and* on the
   *    growth capacity left over by @c i_go. Both residuals are zero as soon
   *    as the cell overflows or ferments, so acetate is never produced and
   *    consumed in the same step.
   *



   *
   * As in SimpleAcetate the internal currency is the elongation rate
   * @c a [m/s], converted to a mass flux through the linear density.
   */
  struct Acetate
  {
    using uniform_weight = std::true_type;
    using Self = Acetate;
    using FloatType = float;
    using Config = std::nullopt_t;

    enum class particle_var : int // NOLINT
    {
      length = 0,
      l_max,
      a_p,
      a_max,
      // Export
      a_e,
      a_go,
      a_fe,
      a_ao,
      // Contributions
      phi_s,
      phi_o2,
      phi_a,
      __COUNT__
    };

    static constexpr std::size_t n_var
        = INDEX_FROM_ENUM(particle_var::__COUNT__);

    static constexpr std::string_view name = "acetate";
    using SelfParticle = MC::ParticlesModel<Self::n_var, Self::FloatType>;

    /// Liquid species: glucose, dioxygen, acetate
    static constexpr std::size_t n_c = 3;

    using SelfContribs = MC::ParticlesContribs<Self::n_c, Self::FloatType>;

    /**
     * @note Plain enumerators, not MODEL_CONSTANT: they are used as runtime
     * subscripts of Kokkos::Array in the device kernels and nvcc rejects a
     * `static constexpr` scalar there ("undefined in device code"). An
     * enumerator has no storage, so it is always usable on the device.
     */
    /// Species indices, O2 has to stay at index 1 (gas/liquid transfer)
    enum species_index : std::size_t
    {
      S = 0,
      O2 = 1,
      Ac = 2,
      N_N = n_c
    };

    /// Reaction indices
    enum reaction_index : std::size_t
    {
      i_go = 0, // glucose oxidation
      i_fe = 1, // glucose fermentation
      i_ov = 2, // acetate overflow
      i_ao = 3, // acetate oxidation
      N_R = 4
    };

    using reaction_rates = Kokkos::Array<FloatType, N_R>;
    using uptake_rates = Kokkos::Array<FloatType, N_N>;

    MODEL_CONSTANT FloatType a_max_m = 2e-6 / 3600.; // m/s
    MODEL_CONSTANT FloatType l_max_m = 2e-6;         // m
    MODEL_CONSTANT FloatType l_min_m = l_max_m / 2.; // m
    MODEL_CONSTANT FloatType d_m = 0.6e-6;           // m
    MODEL_CONSTANT FloatType lin_density
        = c_linear_density(static_cast<FloatType>(1000), d_m); // kg/m

    /// Monod affinity constants [kg/m^3]
    MODEL_CONSTANT FloatType k_s = 1e-3;
    /// Anane et. al 2017 (Biochem. Eng. J). The sharp switch this gives is
    /// what the uptake cap of @ref positivity makes affordable explicitly.
    MODEL_CONSTANT FloatType k_o = 1e-5;
    MODEL_CONSTANT FloatType k_a = 1e-4;

    // Stoichiometry, growth yields

    MODEL_CONSTANT FloatType y_x_go = 0.5F;  ///< gX/gS
    MODEL_CONSTANT FloatType y_x_fe = 0.15F; ///< gX/gS
    MODEL_CONSTANT FloatType y_x_ov = 0.0F;  ///< overflow does not grow
    MODEL_CONSTANT FloatType y_x_ao = 0.33F; ///< gX/gAc

    /// gO2/gS, 0.384
    MODEL_CONSTANT FloatType y_o2_go = static_cast<FloatType>(
        Stoichiometry::y_o2_ox(Stoichiometry::m_cmol_glucose,
                               Stoichiometry::gamma_glucose,
                               y_x_go));
    /// gO2/gAc
    MODEL_CONSTANT FloatType y_o2_ao = static_cast<FloatType>(
        Stoichiometry::y_o2_ox(Stoichiometry::m_cmol_acetate,
                               Stoichiometry::gamma_acetate,
                               y_x_ao));

    /// gAc/gS produced by fermentation, 0.808
    MODEL_CONSTANT FloatType y_ac_fe = static_cast<FloatType>(
        Stoichiometry::y_ac_fer(Stoichiometry::m_cmol_glucose,
                                Stoichiometry::gamma_glucose,
                                y_x_fe));
    /// gAc/gS produced by overflow. C6H12O6 -> 3 C2H4O2 balances as is, so 1
    MODEL_CONSTANT FloatType y_ac_ov = static_cast<FloatType>(
        Stoichiometry::y_ac_fer(Stoichiometry::m_cmol_glucose,
                                Stoichiometry::gamma_glucose,
                                y_x_ov));

    // A reaction cannot need a negative amount of O2 or acetate, nor build
    // more carbon than its substrate carries.
    static_assert(y_o2_go > 0.F
                      && y_x_go <= Stoichiometry::m_cmol_biomass
                                       / Stoichiometry::m_cmol_glucose,
                  "i_go: y_x too high, negative O2 or CO2");
    static_assert(y_o2_ao > 0.F
                      && y_x_ao <= Stoichiometry::m_cmol_biomass
                                       / Stoichiometry::m_cmol_acetate,
                  "i_ao: y_x too high, negative O2 or CO2");
    static_assert(y_ac_fe > 0.F && y_ac_fe < y_ac_ov,
                  "i_fe: y_x too high, cannot grow and excrete as much acetate "
                  "as the pure overflow");

    /// Oxygen uptake capacity, as a multiple of the stoichiometric demand of a
    /// fully respired glucose uptake. Above 1 the cell is never limited by its
    /// own uptake machinery, only by the dissolved oxygen available.
    /// @note the fermentation onset does not depend on @c y_o2_go: the
    /// respirable glucose is min(phi_S, f_o2_max * phi_s_max * monod(O2)), the
    /// yield cancels between the capacity and the phi[O2]/y_o2_go bound.
    MODEL_CONSTANT FloatType f_o2_max = 1.2;
    /// Acetate uptake capacity relative to the glucose one
    MODEL_CONSTANT FloatType f_ac_max = 1. / 3.;

    MODEL_CONSTANT FloatType phi_s_max = a_max_m * lin_density / y_x_go; // kg/s

    /// Highest biomass concentration the reactor is expected to reach
    /// [kg/m^3]. Sizes the cap below, no influence on the metabolism;
    /// overestimating it only makes the cap engage a little earlier.
    MODEL_CONSTANT FloatType x_c_max = 30.;
    /// One explicit step removes at most 1/f_pos of the local concentration
    MODEL_CONSTANT FloatType f_pos = 2.;
    /// Only there so that d_t = 0 does not divide by zero
    MODEL_CONSTANT FloatType d_t_floor = 1e-30;

    /// Lowest admissible growth capacity [m/s]. a_p has to stay *strictly*
    /// positive: it is the only thing that lets a starved lineage grow again,
    /// and log(a_p) is taken at every division.
    MODEL_CONSTANT FloatType a_p_min = 1e-3F * a_max_m;

    /// Regression of the newborn capacity toward its own a_max at division.
    /// a_e <= a_p holds by construction, so a redistribution centred on the
    /// performance of the mother alone (0) is a one way ratchet: the capacity
    /// of a lineage can only decrease and collapses to 0 over the generations.
    MODEL_CONSTANT FloatType a_p_recovery = 0.25F;

    MODEL_CONSTANT auto l_max_dist
        = MC::Distributions::TruncatedNormal<FloatType>(
            l_max_m, l_max_m / 10., l_max_m * 0.7, 1.3 * l_max_m);

    MODEL_CONSTANT auto l_dist = MC::Distributions::TruncatedNormal<FloatType>(
        l_max_m * 0.75, l_max_m / 10., 0.7 * l_min_m, l_max_m * 1.3);

    KOKKOS_INLINE_FUNCTION static void init(const MC::pool_type& random_pool,
                                            std::size_t idx,
                                            const SelfParticle& arr);

    KOKKOS_INLINE_FUNCTION static MC::Status
    update([[maybe_unused]] const MC::pool_type& random_pool,
           FloatType d_t,
           std::size_t idx,
           const SelfParticle& arr,
           const SelfContribs& arr_contribs,
           std::size_t position_index,
           const MC::LocalConcentration& c);

    KOKKOS_INLINE_FUNCTION static void
    division(const MC::pool_type& random_pool,
             std::size_t idx,
             std::size_t idx2,
             const SelfParticle& arr,
             const SelfParticle& buffer_arr);

    /**
     * @brief Distribute the available uptake fluxes over the four reactions.
     * @param phi   uptake potentials [kg/s] indexed by species
     * @param nu_p  potential growth rate of the cell [kgX/s]
     * @param nu    out: growth rate effectively realised [kgX/s]
     * @return the reaction fluxes [kg of substrate/s]
     */
    KOKKOS_INLINE_FUNCTION static reaction_rates
    metabolism(const uptake_rates& phi, FloatType nu_p, FloatType& nu);

    KOKKOS_INLINE_FUNCTION static double
    mass(std::size_t idx, const SelfParticle& arr)
    {
      return GET_PROPERTY(Self::particle_var::length) * lin_density;
    }

    static std::array<std::string_view, Self::n_var>
    names()
    {
      return {
        "length", "l_max", "a_p",   "a_max",  "a_e",   "a_go",
        "a_fe",   "a_ao",  "phi_s", "phi_o2", "phi_a",
      };
    }

    static std::array<std::string_view, Self::n_c>
    species()
    {
      return { "glucose", "O2", "acetate" };
    }
  };

  CHECK_MODEL(Acetate)

  KOKKOS_INLINE_FUNCTION void
  Acetate::init([[maybe_unused]] const MC::pool_type& random_pool,
                std::size_t idx,
                const SelfParticle& arr)
  {

    MODEL_CONSTANT auto a_max_dist
        = MC::Distributions::TruncatedNormal<FloatType>(
            a_max_m, a_max_m / 2., 0.5 * a_max_m, a_max_m * 1.5);

    static constexpr auto ld = l_dist;
    static constexpr auto lm = l_dist;
    auto gen = random_pool.get_state();
    GET_PROPERTY(particle_var::length) = ld.draw(gen);
    GET_PROPERTY(particle_var::l_max) = lm.draw(gen);
    GET_PROPERTY(particle_var::a_p) = a_max_m / 2.;
    GET_PROPERTY(particle_var::a_max) = a_max_dist.mean();
    random_pool.free_state(gen);

    GET_PROPERTY(particle_var::a_e) = 0.;
    GET_PROPERTY(particle_var::a_go) = 0.;
    GET_PROPERTY(particle_var::a_fe) = 0.;
    GET_PROPERTY(particle_var::a_ao) = 0.;
    GET_PROPERTY(particle_var::phi_s) = 0.;
    GET_PROPERTY(particle_var::phi_o2) = 0.;
    GET_PROPERTY(particle_var::phi_a) = 0.;
  }

  KOKKOS_INLINE_FUNCTION Acetate::reaction_rates
  Acetate::metabolism(const uptake_rates& phi,
                      const FloatType nu_p,
                      FloatType& nu)
  {
    constexpr auto zero = static_cast<FloatType>(0);
    reaction_rates r = {}; // zero initialised

    // i_go: respiration, highest priority. Bounded by the glucose, by the
    // oxygen and by the
    // growth capacity.
    const FloatType s_ox_max = Kokkos::min(phi[S], phi[O2] / y_o2_go); // kgS/s
    r[i_go] = Kokkos::min(s_ox_max, nu_p / y_x_go);
    nu = r[i_go] * y_x_go;

    // i_ov: overflow. Glucose that oxygen could have supported but that the
    // growth capacity cannot handle, excreted as acetate. a_p driven branch.
    r[i_ov] = Kokkos::max(s_ox_max - r[i_go], zero);

    // i_ao: acetate re-uptake, on the oxygen and on the growth capacity left
    // over by i_go.  so
    // acetate is never produced and consumed at the same time.
    const FloatType o2_res = Kokkos::max(phi[O2] - (r[i_go] * y_o2_go), zero);
    r[i_ao] = Kokkos::min(Kokkos::min(phi[Ac], o2_res / y_o2_ao),
                          Kokkos::max(nu_p - nu, zero) / y_x_ao);
    nu += r[i_ao] * y_x_ao;

    // i_fe: fermentation, lowest priority.
    const FloatType s_res = Kokkos::max(phi[S] - r[i_go] - r[i_ov], zero);
    r[i_fe] = Kokkos::min(s_res, Kokkos::max(nu_p - nu, zero) / y_x_fe);
    nu += r[i_fe] * y_x_fe;

    KOKKOS_ASSERT(r[i_go] >= zero);
    KOKKOS_ASSERT(r[i_fe] >= zero);
    KOKKOS_ASSERT(r[i_ov] >= zero);
    KOKKOS_ASSERT(r[i_ao] >= zero);
    KOKKOS_ASSERT(nu >= zero && nu <= nu_p * static_cast<FloatType>(1.001));
    // No reaction may draw more than the uptake made available: this is what
    // carries the uptake cap over to the liquid contributions
    KOKKOS_ASSERT(r[i_go] + r[i_fe] + r[i_ov]
                  <= phi[S] * static_cast<FloatType>(1.001));
    KOKKOS_ASSERT((r[i_go] * y_o2_go) + (r[i_ao] * y_o2_ao)
                  <= phi[O2] * static_cast<FloatType>(1.001));
    KOKKOS_ASSERT(r[i_ao] <= phi[Ac] * static_cast<FloatType>(1.001));

    return r;
  }

  KOKKOS_INLINE_FUNCTION MC::Status
  Acetate::update([[maybe_unused]] const MC::pool_type& random_pool,
                  FloatType d_t,
                  std::size_t idx,
                  const SelfParticle& arr,
                  const SelfContribs& arr_contribs,
                  const std::size_t position_index,
                  const MC::LocalConcentration& c)
  {
    const auto c_s = GET_CLAMPED_CONCENTRATION_CAST(FloatType, S);
    const auto c_o2 = GET_CLAMPED_CONCENTRATION_CAST(FloatType, O2);
    const auto c_ac = GET_CLAMPED_CONCENTRATION_CAST(FloatType, Ac);

    const FloatType phi_s_max_p
        = GET_PROPERTY(particle_var::a_max) * lin_density / y_x_go;
    const FloatType phi_o2_max_p = f_o2_max * y_o2_go * phi_s_max_p;
    const FloatType phi_ac_max_p = f_ac_max * phi_s_max_p;

    // Largest uptake one cell may take from the compartment during this step
    // without the explicit liquid update going negative
    const FloatType share = (GET_PROPERTY(particle_var::length) * lin_density)
                            / ((f_pos * x_c_max * d_t) + d_t_floor);

    // Monod uptake potentials [kg/s].
    uptake_rates phi{};
    phi[S] = phi_s_max_p * c_s / (c_s + k_s);
    phi[O2] = Kokkos::min(phi_o2_max_p * c_o2 / (c_o2 + k_o), share * c_o2);
    phi[Ac] = Kokkos::min(phi_ac_max_p * c_ac / (c_ac + k_a), share * c_ac);

    // Growth capacity of the individual, inherited at division and kept in
    // ]0, a_max] by division()
    const FloatType nu_p = GET_PROPERTY(particle_var::a_p) * lin_density;

    FloatType nu = 0.F;
    const reaction_rates r = metabolism(phi, nu_p, nu);

    const FloatType inv_lin_density = static_cast<FloatType>(1) / lin_density;
    const FloatType a_e = nu * inv_lin_density;

    GET_PROPERTY(particle_var::a_e) = a_e;
    GET_PROPERTY(particle_var::a_go) = r[i_go] * y_x_go * inv_lin_density;
    GET_PROPERTY(particle_var::a_fe) = r[i_fe] * y_x_fe * inv_lin_density;
    GET_PROPERTY(particle_var::a_ao) = r[i_ao] * y_x_ao * inv_lin_density;

    GET_PROPERTY(particle_var::length) += d_t * a_e;

    // Glucose: consumed by respiration, fermentation and overflow
    const FloatType phi_s = -(r[i_go] + r[i_fe] + r[i_ov]);
    // Oxygen: only consumed, it is fed to the liquid by the gas/liquid transfer
    const FloatType phi_o2 = -((r[i_go] * y_o2_go) + (r[i_ao] * y_o2_ao));
    // Acetate: produced by fermentation (O2 limited) and by overflow (a_p
    // limited), consumed by acetate oxidation
    const FloatType phi_a
        = ((r[i_fe] * y_ac_fe) + (r[i_ov] * y_ac_ov)) - r[i_ao];

    GET_PROPERTY(particle_var::phi_s) = phi_s;
    GET_PROPERTY(particle_var::phi_o2) = phi_o2;
    GET_PROPERTY(particle_var::phi_a) = phi_a;

    GET_CONTRIBS(S) = phi_s;
    GET_CONTRIBS(O2) = phi_o2;
    GET_CONTRIBS(Ac) = phi_a;

    return check_div(GET_PROPERTY(Self::particle_var::length),
                     GET_PROPERTY(Self::particle_var::l_max));
  }

  KOKKOS_INLINE_FUNCTION void
  Acetate::division([[maybe_unused]] const MC::pool_type& random_pool,
                    std::size_t idx,
                    std::size_t idx2,
                    const SelfParticle& arr,
                    const SelfParticle& buffer_arr)
  {
    Kokkos::View<FloatType**,
                 ComputeSpace::array_layout,
                 Kokkos::MemoryTraits<Kokkos::MemoryTraitsFlags::Restrict>>
        buffer_e = buffer_arr;

    const FloatType current_l = GET_PROPERTY(particle_var::length);
    const FloatType new_current_length = current_l / 2.F;
    GET_PROPERTY(particle_var::length) = new_current_length;
    constexpr auto binf = INDEX_FROM_ENUM(particle_var::length);
    constexpr auto bsup = INDEX_FROM_ENUM(particle_var::a_e);
    for (auto i = binf; i < bsup; ++i)
    {
      COPY_PROPERTY_TO(i, idx2, buffer_e);
    }

    // Local copies: nvcc rejects a static constexpr member passed by reference
    // to Kokkos::min/max/clamp ("undefined in device code")
    const FloatType ap_min = a_p_min;
    const FloatType a_max_i
        = Kokkos::max(GET_PROPERTY(particle_var::a_max), ap_min);

    // Effective elongation of the mother, floored: a fully starved cell would
    // otherwise give log(0)
    const FloatType a_e_i
        = Kokkos::clamp(GET_PROPERTY(particle_var::a_e), ap_min, a_max_i);

    // Centre of the redistribution. a_e <= a_p always holds, so centring on
    // a_e alone makes the capacity of a lineage a decreasing sequence that
    // converges to 0; the pull toward a_max keeps a_p away from 0 and lets a
    // starved lineage recover when the environment improves.
    const FloatType a_target = a_e_i + (a_p_recovery * (a_max_i - a_e_i));

    auto gen = random_pool.get_state();
    const double sigma = 0.2;
    const double average
        = Kokkos::log(static_cast<double>(a_target)) - sigma * sigma / 2;
    const auto dist = MC::Distributions::LogNormal<double>(average, sigma);

    // Clamped to [a_p_min, a_max]: a null capacity freezes the sub-lineage for
    // the rest of the run
    const auto gen1 = Kokkos::clamp(
        static_cast<FloatType>(dist.draw(gen)), ap_min, a_max_i);
    const auto gen2 = Kokkos::clamp(
        static_cast<FloatType>(dist.draw(gen)), ap_min, a_max_i);

    static constexpr auto local_l = l_max_dist;
    const FloatType lmax1 = local_l.draw(gen);
    const FloatType lmax2 = local_l.draw(gen);

    GET_PROPERTY(particle_var::a_p) = gen1;
    GET_PROPERTY(particle_var::l_max) = lmax1;
    GET_PROPERTY_FROM(idx2, buffer_arr, particle_var::a_p) = gen2;
    GET_PROPERTY_FROM(idx2, buffer_arr, particle_var::l_max) = lmax2;

    random_pool.free_state(gen);
  }

} // namespace Models

#endif
