#ifndef __MODEL_MONOD_HPP__
#define __MODEL_MONOD_HPP__

#include "mc/macros.hpp"
#include "models/utils.hpp"
#include <Kokkos_MathematicalConstants.hpp>
#include <mc/alias.hpp>
#include <mc/prng/prng_extension.hpp>
#include <mc/traits.hpp>
#include <optional>

namespace Models
{

  /**
   * @brief Simplified Monod model for glucose consumption and biomass growth
   *
   * This model links glucose consumption and biomass growth using the Monod
   * equation. It models biomass growth through cell elongation for rod-shaped
   * strains, and cell division is assumed to be deterministic. The primary
   * purpose of this model is for validation and to demonstrate how the
   * framework operates without requiring the typical formation. Instead, it
   * focuses on the core dynamics of glucose consumption and biomass growth
   * through the Monod law.
   */
  struct Monod
  {
    /**
     * @brief Enumeration for the Monod model variables
     *
     * This enum represents the different variables used in the Monod model.
     * These variables do not strictly correspond to the state vector (xi), as
     * some are auxiliary variables (e.g., phi_s_c) not part of the primary
     * state. The Monod factor is used to convert the growth rate into
     * elongation rates.
     */
    enum class particle_var : int
    {
      l = 0, ///< Length of the cell
      l_max, ///< Maximum cell length (assumed constant for all cells)
      mu_p,  ///< Maximum specific growth rate (potential growth rate)
      mue,   ///< Effective growth rate (used for export purposes only)
      _init_only_cell_lenghtening, ///< Factor for converting growth rate into
                                   ///< elongation rate
      phi_s_c,                     ///< Instantaneous glucose consumption rate
      __COUNT__ ///< Helper for determining the size of the variable list
    };

    // clang-format off
    static constexpr std::size_t n_var = INDEX_FROM_ENUM(particle_var::__COUNT__);
    static constexpr std::string_view name = "monod";
    using uniform_weight = std::true_type; // Using type alias
    using Self = Monod;
    using FloatType = float;
    using SelfParticle = MC::ParticlesModel<Self::n_var, Self::FloatType>;
    using Config = std::nullopt_t;
    

    /** 
    * @brief Model constants used in biomass growth and cell elongation
    * 
    * These constants define various parameters for the model, including maximum growth rates,
    * cell size, and other metabolic factors. The units are specified in meters (m) unless otherwise indicated.
    */
    MODEL_CONSTANT FloatType y_s_x = 2;                ///< Specific yield of biomass per unit substrate (m)
    MODEL_CONSTANT FloatType mu_max = 0.77 / 3600.;    ///< Maximum specific growth rate (1/s), converted from per hour
    MODEL_CONSTANT FloatType tau_meta = 1. / mu_max;   ///< Metabolic time constant (s), inverse of max growth rate
    MODEL_CONSTANT FloatType l_max_m = 2e-6;           ///< Maximum cell length (m)
    MODEL_CONSTANT FloatType l_min_m = l_max_m / 2.;   ///< Minimum cell length (m), half of maximum length
    MODEL_CONSTANT FloatType k_s = 1e-3;               ///< Monod constant for substrate concentration (m)
    MODEL_CONSTANT FloatType d_m = 0.6e-6;             ///< Cell diameter (m)
    MODEL_CONSTANT FloatType lin_density = 
        c_linear_density(static_cast<FloatType>(1000), d_m); ///< Linear density of the biomass (kg/m), calculated from cell diameter
    // clang-format on

    MODEL_CONSTANT auto initial_length_dist =
        MC::Distributions::TruncatedNormal<FloatType>(
            l_max_m * 0.75, l_max_m * 0.75 / 4, l_min_m, l_max_m);

    KOKKOS_INLINE_FUNCTION static void
    init([[maybe_unused]] const MC::KPRNG::pool_type& random_pool,
         [[maybe_unused]] std::size_t idx,
         [[maybe_unused]] const SelfParticle& arr)
    {
      // Helper cause, initial_length_dist may not be captured in cuda context
      // (FIXME)
      static constexpr auto local_dist = initial_length_dist;

      const auto [l0, mu] = MC::sample_random_variables(
          random_pool,
          [](auto gen)
          {
            const auto l0 = local_dist.draw(
                gen); // Get initial fro given distribution (normal)
            const auto mu = mu_max; // Set maximal for each cell mu_p
            return std::make_tuple(l0, mu);
          });

      GET_PROPERTY(particle_var::l) = l0;
      GET_PROPERTY(particle_var::l_max) =
          l_max_m; // Set the same length for everyone
      GET_PROPERTY(particle_var::mu_p) = mu_max;

      // born at l_max/2, divides at l_max
      // Need to be changed if division criteria is changed
      constexpr auto dl = l_max_m / 2.;

      // Monod factor to convert growth rate into elongation, it assumes
      // doubling in mass during generation time G=ln(2)/mu factor= DL/ln(2)
      // comes from mu=ln(2)/G
      GET_PROPERTY(particle_var::_init_only_cell_lenghtening) =
          dl / Kokkos::numbers::ln2;
    }

    KOKKOS_INLINE_FUNCTION static double
    mass([[maybe_unused]] std::size_t idx,
         [[maybe_unused]] const SelfParticle& arr)
    {
      return GET_PROPERTY(Self::particle_var::l) * lin_density; // Rod-shaped
    }

    KOKKOS_INLINE_FUNCTION static MC::Status
    update([[maybe_unused]] const MC::KPRNG::pool_type& random_pool,
           [[maybe_unused]] FloatType d_t,
           [[maybe_unused]] std::size_t idx,
           [[maybe_unused]] const SelfParticle& arr,
           [[maybe_unused]] const MC::LocalConcentration& c)
    {
      const FloatType s =
          static_cast<FloatType>(Kokkos::max(0., c(0))); // Bounded
      const FloatType mu = mu_max * s / (k_s + s); // Instantaneous mu from
                                                   // Monod

      // Efffective growth rate
      const FloatType mu_eff =
          Kokkos::min(GET_PROPERTY(Self::particle_var::mu_p), mu);

      // Lengthening
      GET_PROPERTY(Self::particle_var::l) +=
          d_t * (mu_eff *
                 GET_PROPERTY(Self::particle_var::_init_only_cell_lenghtening));

      // Growth rate
      GET_PROPERTY(Self::particle_var::mu_p) +=
          d_t * (1.0 / tau_meta) *
          (mu - GET_PROPERTY(Self::particle_var::mu_p));

      // Store only for being exported
      GET_PROPERTY(Self::particle_var::mue) = mu_eff;

      // Contributions
      GET_PROPERTY(Self::particle_var::phi_s_c) =
          -mu_eff * y_s_x * static_cast<FloatType>(mass(idx, arr));

      return check_div(GET_PROPERTY(Self::particle_var::l),
                       GET_PROPERTY(Self::particle_var::l_max));
    }

    KOKKOS_INLINE_FUNCTION static void
    division([[maybe_unused]] const MC::KPRNG::pool_type& random_pool,
             [[maybe_unused]] std::size_t idx,
             [[maybe_unused]] std::size_t idx2,
             [[maybe_unused]] const SelfParticle& arr,
             [[maybe_unused]] const SelfParticle& buffer_arr)
    {

      //_init_only_cell_lenghtening doesnt need to be redistributed bcause
      // division criteria is the same for each cell then dl remains the same

      // Even if division is deterministic, keep generate way to divide

      /** @note  Important to note that new_current_length is not l_max/2 but
      l/2 This is completly different but as growth is slow and time step has to
      be small it's ok */

      const FloatType current_l = GET_PROPERTY(particle_var::l);
      const FloatType new_current_length = current_l / 2.F;
      GET_PROPERTY(particle_var::l) = new_current_length;

      const FloatType mu_p_value = GET_PROPERTY(particle_var::mu_p);
      const FloatType cell_lenghtening_value =
          GET_PROPERTY(particle_var::_init_only_cell_lenghtening);

      GET_PROPERTY_FROM(idx2, buffer_arr, particle_var::l) = new_current_length;
      GET_PROPERTY_FROM(idx2, buffer_arr, particle_var::l_max) = l_max_m;
      GET_PROPERTY_FROM(idx2, buffer_arr, particle_var::mu_p) = mu_p_value;
      GET_PROPERTY_FROM(
          idx2, buffer_arr, particle_var::_init_only_cell_lenghtening) =
          cell_lenghtening_value;
    }

    static MC::ContribIndexBounds get_bounds()
    {
      int begin = INDEX_FROM_ENUM(Self::particle_var::phi_s_c);
      return {.begin = begin, .end = begin + 1};
    }

    static std::vector<std::string_view> names()
    {
      return {

          "length", "mu", "mu_eff"};
    }

    static std::vector<std::size_t> get_number()
    {
      return {INDEX_FROM_ENUM(particle_var::l),
              INDEX_FROM_ENUM(particle_var::mu_p),
              INDEX_FROM_ENUM(particle_var::mue)};
    }
  };

  CHECK_MODEL(Monod)
} // namespace Models
#endif
