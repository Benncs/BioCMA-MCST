
#ifndef __FIXED_LENGTH_MODEL_HPP__
#define __FIXED_LENGTH_MODEL_HPP__

#include "common/traits.hpp"
#include "mc/macros.hpp"
#include "models/utils.hpp"
#include <mc/prng/prng_extension.hpp>
#include <mc/traits.hpp>
#include <optional>
#include <string_view>

namespace Models
{
  template <FloatingPointType F>
  static F consteval _get_phi_s_max(F density, F dl)
  {
    // dl and density must be same unit, dl*density -> mass and y is mass yield
    return (dl * density) * 0.5;
  }
  struct FixedLength
  {
    using uniform_weight = std::true_type; // Using type alias
    using Self = FixedLength;
    using FloatType = float;
    using Config = std::nullopt_t;

    enum class particle_var : int
    {
      length = 0,
      l_max,
      phi_s,
      COUNT
    };

    static constexpr std::size_t n_var = INDEX_FROM_ENUM(particle_var::COUNT);
    static constexpr std::string_view name = "fixed-length";
    using SelfParticle = MC::ParticlesModel<Self::n_var, Self::FloatType>;

    MODEL_CONSTANT FloatType l_dot_max = 4e-6 / 3600.; // m
    MODEL_CONSTANT FloatType l_max_m = 2e-6;           // m
    MODEL_CONSTANT FloatType k = 1e-3;                 // m
    MODEL_CONSTANT FloatType d_m = 0.6e-6;             // m
    MODEL_CONSTANT FloatType l_min_m = l_max_m / 2.;   // m
    MODEL_CONSTANT FloatType lin_density =
        c_linear_density(static_cast<FloatType>(1000), d_m);

    MODEL_CONSTANT FloatType phi_s_max =
        _get_phi_s_max<FloatType>(lin_density, l_dot_max * 1.1); // kgS/s

    KOKKOS_INLINE_FUNCTION static void
    init(const MC::KPRNG::pool_type& random_pool,
         std::size_t idx,
         const SelfParticle& arr);

    KOKKOS_INLINE_FUNCTION static MC::Status
    update(const MC::KPRNG::pool_type& random_pool,
           FloatType d_t,
           std::size_t idx,
           const SelfParticle& arr,
           const MC::LocalConcentration& c);

    KOKKOS_INLINE_FUNCTION static void
    division(const MC::KPRNG::pool_type& random_pool,
             std::size_t idx,
             std::size_t idx2,
             const SelfParticle& arr,
             const SelfParticle& buffer_arr);

    KOKKOS_INLINE_FUNCTION static void
    contribution(std::size_t idx,
                 std::size_t position,
                 double weight,
                 const SelfParticle& arr,
                 const MC::ContributionView& contributions);

    KOKKOS_INLINE_FUNCTION static double mass(std::size_t idx,
                                              const SelfParticle& arr)
    {
      return GET_PROPERTY(Self::particle_var::length) * lin_density;
    }
    static std::vector<std::string_view> names()
    {
      return {

          "length",
      };
    }
    static std::vector<std::size_t> get_number()
    {
      return {INDEX_FROM_ENUM(particle_var::length)};
    }
  };

  CHECK_MODEL(FixedLength)

  KOKKOS_INLINE_FUNCTION void
  FixedLength::init([[maybe_unused]] const MC::KPRNG::pool_type& random_pool,
                    std::size_t idx,
                    const SelfParticle& arr)
  {

    constexpr auto l_initial_dist = MC::Distributions::TruncatedNormal<float>(
        l_min_m, l_min_m / 7., 0., static_cast<double>(1));

    auto gen = random_pool.get_state();
    auto linit = l_initial_dist.draw(
        gen); // l_min_m; // gen.drand(l_min_m * .8, l_max_m);
    GET_PROPERTY(particle_var::l_max) =
        l_max_m; // gen.drand(l_max_m * 0.8, l_max_m * 1.2);
    random_pool.free_state(gen);
    GET_PROPERTY(particle_var::length) = linit;
    GET_PROPERTY(particle_var::phi_s) = 0.;
  }

  KOKKOS_INLINE_FUNCTION MC::Status
  FixedLength::update([[maybe_unused]] const MC::KPRNG::pool_type& random_pool,
                      FloatType d_t,
                      std::size_t idx,
                      const SelfParticle& arr,
                      const MC::LocalConcentration& c)
  {
    const FloatType g = c(0) / (k + c(0));
    FloatType phi_s = phi_s_max * g;
    const auto ldot = l_dot_max * g;

    GET_PROPERTY(FixedLength::particle_var::phi_s) = phi_s * 2.;
    GET_PROPERTY(FixedLength::particle_var::length) += d_t * ldot;

    return (GET_PROPERTY(Self::particle_var::length) >=
            GET_PROPERTY(particle_var::l_max))
               ? MC::Status::Division
               : MC::Status::Idle;
  }

  KOKKOS_INLINE_FUNCTION void
  FixedLength::division(const MC::KPRNG::pool_type& random_pool,
                        std::size_t idx,
                        std::size_t idx2,
                        const SelfParticle& arr,
                        const SelfParticle& buffer_arr)
  {
    const FloatType new_current_length =
        GET_PROPERTY(particle_var::length) / 2.F;
    GET_PROPERTY_FROM(idx2, buffer_arr, particle_var::length) =
        new_current_length;
    GET_PROPERTY(particle_var::length) = new_current_length;
    auto gen = random_pool.get_state();
    GET_PROPERTY_FROM(idx2, buffer_arr, particle_var::l_max) =
        gen.drand(l_max_m * 0.8, l_max_m * 1.2);
    random_pool.free_state(gen);
  }

  KOKKOS_INLINE_FUNCTION void
  FixedLength::contribution([[maybe_unused]] std::size_t idx,
                            std::size_t position,
                            double weight,
                            [[maybe_unused]] const SelfParticle& arr,
                            const MC::ContributionView& contributions)
  {
    auto access = contributions.access();
    access(position, 0) +=
        -weight * GET_PROPERTY(FixedLength::particle_var::phi_s); // NOLINT
  }

} // namespace Models

#endif
