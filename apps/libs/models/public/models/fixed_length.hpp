
#ifndef __FIXED_LENGTH_MODEL_HPP__
#define __FIXED_LENGTH_MODEL_HPP__

#include "Kokkos_Clamp.hpp"
#include "Kokkos_Core_fwd.hpp"
#include "Kokkos_Macros.hpp"
#include "common/common.hpp"
#include "common/traits.hpp"
#include "mc/macros.hpp"
#include "models/utils.hpp"
#include <mc/prng/prng_extension.hpp>
#include <mc/traits.hpp>
#include <string_view>

namespace Models
{
  template <FloatingPointType F>
  static F consteval _get_phi_s_max(F density,
                                    F dl,
                                    F glucose_to_biomass_yield = 0.5)
  {
    // dl and density must be same unit, dl*density -> mass and y is mass yield
    return (dl * density) / glucose_to_biomass_yield;
  }
  struct FixedLength
  {
    using uniform_weight = std::true_type;
    using Self = FixedLength;
    using FloatType = float;

    using Config = Kokkos::View<const FloatType*, ComputeSpace>;

    enum class particle_var : int // NOLINT
    {
      length = 0,
      l_max,
      phi_s,
      COUNT
    };

    static constexpr std::size_t n_var = INDEX_FROM_ENUM(particle_var::COUNT);
    static constexpr std::string_view name = "fixed-length";
    using SelfParticle = MC::ParticlesModel<Self::n_var, Self::FloatType>;

    MODEL_CONSTANT FloatType l_dot_max = 2e-6 / 3600.; // m
    MODEL_CONSTANT FloatType l_max_m = 2e-6;           // m
    MODEL_CONSTANT FloatType l_min_m = l_max_m / 2.;   // m
    MODEL_CONSTANT FloatType k = 1e-3;                 // m
    MODEL_CONSTANT FloatType d_m = 0.6e-6;             // m
    MODEL_CONSTANT FloatType lin_density =
        c_linear_density(static_cast<FloatType>(1000), d_m);

    MODEL_CONSTANT FloatType phi_s_max =
        _get_phi_s_max<FloatType>(lin_density, l_dot_max); // kgS/s

    MC::ContribIndexBounds static get_bounds()
    {
      int begin = INDEX_FROM_ENUM(Self::particle_var::phi_s);
      return {.begin = begin, .end = begin + 1};
    }

    static Self::Config get_config(std::size_t n);

    KOKKOS_INLINE_FUNCTION static void
    init(const MC::KPRNG::pool_type& random_pool,
         std::size_t idx,
         const SelfParticle& arr,
         const Config& params);

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
                    const SelfParticle& arr,
                    const Config& config)
  {

    // normal1
    constexpr auto l_initial_dist = MC::Distributions::TruncatedNormal<float>(
        l_max_m * 0.75, l_min_m / 5., l_min_m, l_max_m);

    // normal2 lmin/2 lmin/2
    // normal3 l_min*0.7
    // norma4l
    // constexpr auto l_initial_dist =
    // MC::Distributions::TruncatedNormal<float>(
    //    l_min_m, l_min_m / 2., 0., static_cast<double>(1));

    // Normal1
    // constexpr FloatType lambda =
    //     0.69314718056 / l_max_m;
    //
    //     // lambda=-h/a WARN!!! use abs(lambda)
    //  for exp(-lambda x) sampling
    //     constexpr auto l_initial_dist =

    auto gen = random_pool.get_state();
    // auto linit = l_initial_dist.draw(gen);
    //  FloatType linit = 0.;
    //  do
    //  {
    //    linit = lm / 2. + config().dist.draw(gen);
    //  } while (linit > lm);
    //
    //
    //
    //
    auto linit = config(idx); // linit = Kokkos::min(linit, lm);

    // linit = Kokkos::max(linit, FloatType(lm / 2.));

    //   auto linit = l_initial_dist.draw(
    //       gen); // l_min_m; // gen.drand(l_min_m * .8, l_max_m);
    //
    random_pool.free_state(gen);
    GET_PROPERTY(particle_var::length) = linit;
    // Stochastic 1
    //    auto gen = random_pool.get_state();
    //    GET_PROPERTY(particle_var::l_max) = gen.drand(l_max_m * 0.9, l_max_m
    //    * 1.1); random_pool.free_state(gen);

    // Stochastic 2
    // auto gen = random_pool.get_state();
    // GET_PROPERTY(particle_var::l_max) = gen.drand(l_max_m * 0.8, l_max_m
    // * 1.2); random_pool.free_state(gen);
    GET_PROPERTY(particle_var::l_max) = l_max_m;

    //   GET_PROPERTY(particle_var::length) = l_max_m / 2.;

    GET_PROPERTY(particle_var::phi_s) = 0.;
  }

  KOKKOS_INLINE_FUNCTION MC::Status
  FixedLength::update([[maybe_unused]] const MC::KPRNG::pool_type& random_pool,
                      FloatType d_t,
                      std::size_t idx,
                      const SelfParticle& arr,
                      const MC::LocalConcentration& c)
  {
    const auto s = static_cast<FloatType>(c(0));
    const FloatType g = s / (k + s);
    const FloatType phi_s = phi_s_max * g;
    const FloatType ldot = l_dot_max * g;
    GET_PROPERTY(Self::particle_var::length) += d_t * ldot;
    GET_PROPERTY(Self::particle_var::phi_s) = -phi_s;
    return check_div(GET_PROPERTY(Self::particle_var::length),
                     GET_PROPERTY(Self::particle_var::l_max));
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

    GET_PROPERTY(particle_var::length) = new_current_length;
    GET_PROPERTY(particle_var::l_max) = l_max_m;

    GET_PROPERTY_FROM(idx2, buffer_arr, particle_var::length) =
        new_current_length;
    GET_PROPERTY_FROM(idx2, buffer_arr, particle_var::l_max) = l_max_m;

    // Stochastic 2
    //    auto gen = random_pool.get_state();
    //    GET_PROPERTY(particle_var::l_max) = gen.drand(l_max_m * 0.8, l_max_m
    //    * 1.2); GET_PROPERTY_FROM(idx2, buffer_arr, particle_var::l_max) =
    //        gen.drand(l_max_m * 0.8, l_max_m * 1.2);
    //    random_pool.free_state(gen);

    // Stochastic 1
    // auto gen = random_pool.get_state();
    // GET_PROPERTY(particle_var::l_max) = gen.drand(l_max_m * 0.9, l_max_m
    // * 1.1); GET_PROPERTY_FROM(idx2, buffer_arr, particle_var::l_max) =
    // gen.drand(l_max_m * 0.9, l_max_m * 1.1);
    // random_pool.free_state(gen);
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
