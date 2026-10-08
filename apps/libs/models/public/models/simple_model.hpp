#ifndef __SIMPLE_MODEL_HPP__
#define __SIMPLE_MODEL_HPP__

#include <Kokkos_Core_fwd.hpp>
#include <Kokkos_Macros.hpp>
#include <common/common.hpp>
#include <common/traits.hpp>
#include <mc/alias.hpp>
#include <mc/macros.hpp>
#include <mc/prng/prng_extension.hpp>
#include <mc/traits.hpp>
#include <models/utils.hpp>
#include <optional>
#include <string_view>

namespace Models
{
  /**
   * @brief Rod-shaped cell growing on glucose, dividing on a timer
   *
   * Elongation follows a Monod law on glucose. Division is not triggered by
   * length but by an age threshold drawn per cell, so the population spreads
   * even with a deterministic growth rate.
   */
  struct SimpleModel
  {
    using uniform_weight = std::true_type;
    using Self = SimpleModel;
    using FloatType = float;
    using Config = std::nullopt_t;

    enum class particle_var : int // NOLINT
    {
      length = 0,
      age,
      t_div,
      phi_s,
      __COUNT__
    };

    MODEL_CONSTANT std::size_t n_var = INDEX_FROM_ENUM(particle_var::__COUNT__);

    MODEL_CONSTANT std::size_t n_c = 1;

    MODEL_CONSTANT std::string_view name = "simple";
    using SelfParticle = MC::ParticlesModel<Self::n_var, Self::FloatType>;
    using SelfContribs = MC::ParticlesContribs<Self::n_c, Self::FloatType>;

    MODEL_CONSTANT FloatType l_dot_max = 8 * 2e-10; // m/s
    MODEL_CONSTANT FloatType l_max_m = 5e-6;        // m
    MODEL_CONSTANT FloatType l_min_m = 0.9e-6;      // m
    MODEL_CONSTANT FloatType d_m = 0.6e-6;          // m
    MODEL_CONSTANT FloatType k_s = 1e-3;            // kg/m3
    MODEL_CONSTANT FloatType lin_density
        = c_linear_density(static_cast<FloatType>(1000), d_m);

    MODEL_CONSTANT FloatType phi_s_max
        = _get_phi_s_max<FloatType>(lin_density, l_dot_max); // kgS/s

    MODEL_CONSTANT auto l_dist = MC::Distributions::TruncatedNormal<FloatType>(
        l_min_m, l_min_m / 5., l_min_m * 0.5, l_max_m);

    MODEL_CONSTANT auto t_div_dist
        = MC::Distributions::TruncatedNormal<FloatType>(500.,
                                                        500. / 2.,
                                                        10.,
                                                        1200.); // s

    KOKKOS_INLINE_FUNCTION static void init(const MC::pool_type& random_pool,
                                            std::size_t idx,
                                            const SelfParticle& arr);

    KOKKOS_INLINE_FUNCTION static MC::Status
    update(const MC::pool_type& random_pool,
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

    KOKKOS_INLINE_FUNCTION static double
    mass(std::size_t idx, const SelfParticle& arr)
    {
      return GET_PROPERTY(Self::particle_var::length) * lin_density;
    }

    // age is not exported is it export by default
    static std::vector<std::string_view>
    names()
    {
      return { "length", "t_div", "phi_s" };
    }

    static std::vector<std::size_t>
    get_number()
    {
      return { INDEX_FROM_ENUM(particle_var::length),
               INDEX_FROM_ENUM(particle_var::t_div),
               INDEX_FROM_ENUM(particle_var::phi_s) };
    }

    static std::vector<std::string_view>
    species()
    {
      return { "S" };
    }
  };

  CHECK_MODEL(SimpleModel)

  KOKKOS_INLINE_FUNCTION void
  SimpleModel::init(const MC::pool_type& random_pool,
                    std::size_t idx,
                    const SelfParticle& arr)
  {
    // Local copies, struct scope distributions may not be captured in cuda
    static constexpr auto local_l = l_dist;
    static constexpr auto local_t = t_div_dist;

    auto gen = random_pool.get_state();
    GET_PROPERTY(particle_var::length) = local_l.draw(gen);
    GET_PROPERTY(particle_var::t_div) = local_t.draw(gen);
    random_pool.free_state(gen);

    GET_PROPERTY(particle_var::age) = 0.;
    GET_PROPERTY(particle_var::phi_s) = 0.;
  }

  KOKKOS_INLINE_FUNCTION MC::Status
  SimpleModel::update([[maybe_unused]] const MC::pool_type& random_pool,
                      FloatType d_t,
                      std::size_t idx,
                      const SelfParticle& arr,
                      const SelfContribs& arr_contribs,
                      const std::size_t position_index,
                      const MC::LocalConcentration& c)
  {
    const auto s = GET_CLAMPED_CONCENTRATION_CAST(FloatType, 0);

    const FloatType g = s / (k_s + s);
    const FloatType phi_s = phi_s_max * g;

    GET_PROPERTY(particle_var::length) += d_t * l_dot_max * g;
    GET_PROPERTY(particle_var::age) += d_t;
    GET_PROPERTY(particle_var::phi_s) = phi_s;

    GET_CONTRIBS(0) = -phi_s;

    return (GET_PROPERTY(Self::particle_var::age)
            >= GET_PROPERTY(Self::particle_var::t_div))
               ? MC::Status::Division
               : MC::Status::Idle;
  }

  KOKKOS_INLINE_FUNCTION void
  SimpleModel::division(const MC::pool_type& random_pool,
                        std::size_t idx,
                        std::size_t idx2,
                        const SelfParticle& arr,
                        const SelfParticle& buffer_arr)
  {
    static constexpr auto local_t = t_div_dist;

    const FloatType new_current_length
        = GET_PROPERTY(particle_var::length) / 2.F;

    GET_PROPERTY(particle_var::length) = new_current_length;
    GET_PROPERTY(particle_var::age) = 0.;
    GET_PROPERTY_FROM(idx2, buffer_arr, particle_var::length)
        = new_current_length;
    GET_PROPERTY_FROM(idx2, buffer_arr, particle_var::age) = 0.;

    // Redraw both, without it the whole lineage keeps the same period
    auto gen = random_pool.get_state();
    GET_PROPERTY(particle_var::t_div) = local_t.draw(gen);
    GET_PROPERTY_FROM(idx2, buffer_arr, particle_var::t_div)
        = local_t.draw(gen);
    random_pool.free_state(gen);

    GET_PROPERTY_FROM(idx2, buffer_arr, particle_var::phi_s)
        = GET_PROPERTY(particle_var::phi_s);
  }

} // namespace Models

#endif
