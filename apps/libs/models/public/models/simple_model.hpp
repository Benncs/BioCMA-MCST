#ifndef __SIMPLE_MODEL_HPP__
#define __SIMPLE_MODEL_HPP__

#include "common/traits.hpp"
#include "models/uptake_dyn.hpp"
#include "models/utils.hpp"
#include <mc/prng/prng_extension.hpp>
#include <mc/traits.hpp>
#include <models/uptake.hpp>
#include <optional>
#include <string_view>

namespace Models
{
  template <FloatingPointType F> static F consteval get_phi_s_max(F density, F dl)
  {
    // dl and density must be same unit, dl*density -> mass and y is mass yield
    return (dl * density) * 0.5;
  }
  struct SimpleModel
  {
    using uniform_weight = std::true_type; // Using type alias
    using Self = SimpleModel;
    using FloatType = float;
    using Config = std::nullopt_t;

    enum class particle_var : int
    {
      length = INDEX_FROM_ENUM(Uptakeparticle_var::COUNT),
      age,
      phi_s,
      t_div,
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

    MODEL_CONSTANT FloatType phi_s_max = get_phi_s_max<FloatType>(lin_density, 8 * 2e-10); // kgS/s

    MODEL_CONSTANT FloatType phi_perm_max = phi_s_max / 40.; // kgS/

    MODEL_CONSTANT FloatType frequency_division = 1. / 1000; // Hz

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
      return GET_PROPERTY(SimpleModel::particle_var::length) * lin_density;
    }

    // inline constexpr static std::array<std::string_view, n_var> names()
    // {
    //   constexpr std::size_t ln_var = n_var - Uptake<SimpleModel>::n_var;
    //   constexpr auto _names = concat_arrays<Uptake<SimpleModel>::n_var, ln_var>(
    //       Uptake<SimpleModel>::names(), {"length", "age", "phi_s", "t_div"});

    //   return _names;
    // }


    
    
  };

  CHECK_MODEL(SimpleModel)

  KOKKOS_INLINE_FUNCTION void
  SimpleModel::init([[maybe_unused]] const MC::KPRNG::pool_type& random_pool,
                    std::size_t idx,
                    const SelfParticle& arr)
  {
    MODEL_CONSTANT auto division_time_d = MC::Distributions::TruncatedNormal<FloatType>(
        500, 500./2.,10,1200); //-Kokkos::log(random_number) / frequency_division;
    MODEL_CONSTANT auto l_distribution = MC::Distributions::TruncatedNormal<FloatType>(
        l_min_m, l_min_m / 5., l_min_m * 0.5, l_max_m);
    auto gen = random_pool.get_state();
    arr(idx, static_cast<int>(particle_var::length)) = l_distribution.draw(gen);
    GET_PROPERTY(SimpleModel::particle_var::t_div) = division_time_d.draw(gen);
    random_pool.free_state(gen);
    arr(idx, static_cast<int>(particle_var::age)) = 0;
    arr(idx, static_cast<int>(particle_var::phi_s)) = 0;

    Uptake<UptakeDefault<typename Self::FloatType>,SimpleModel>::init(random_pool, idx, arr);
  }

  KOKKOS_INLINE_FUNCTION MC::Status SimpleModel::update([[maybe_unused]]const MC::KPRNG::pool_type& random_pool,
                                                        FloatType d_t,
                                                        std::size_t idx,
                                                        const SelfParticle& arr,
                                                        const MC::LocalConcentration& c)
  {

    const auto phi_s = Uptake<UptakeDefault<typename Self::FloatType>,SimpleModel>::uptake_step(phi_s_max, d_t, idx, arr, c);

    GET_PROPERTY(SimpleModel::particle_var::age) += d_t;
    GET_PROPERTY(SimpleModel::particle_var::phi_s) = phi_s;

    return (GET_PROPERTY(SimpleModel::particle_var::t_div) <=
            GET_PROPERTY(SimpleModel::particle_var::age))
               ? MC::Status::Division
               : MC::Status::Idle;
  }

  KOKKOS_INLINE_FUNCTION void SimpleModel::division(const MC::KPRNG::pool_type& random_pool,
                                                    std::size_t idx,
                                                    std::size_t idx2,
                                                    const SelfParticle& arr,
                                                    const SelfParticle& buffer_arr)
  {
    // for(size_t i=0;i<n_var;++i)
    // {
    //   buffer_arr(idx,i)=arr(idx,i); //basic copy
    // }
    MODEL_CONSTANT auto division_time_d = MC::Distributions::TruncatedNormal<FloatType>(
        500., 500./2.,10,1200); //-Kokkos::log(random_number) / frequency_division;

    const FloatType new_current_length =
        arr(idx, static_cast<int>(particle_var::length)) / static_cast<FloatType>(2.);
    buffer_arr(idx2, static_cast<int>(particle_var::length)) = new_current_length;
    arr(idx, static_cast<int>(particle_var::length)) = new_current_length;

    arr(idx, static_cast<int>(particle_var::age)) = 0;
    buffer_arr(idx2, static_cast<int>(particle_var::age)) = 0;

    auto gen = random_pool.get_state();
    GET_PROPERTY_FROM(idx2, buffer_arr, SimpleModel::particle_var::t_div) =
        division_time_d.draw(gen);
    GET_PROPERTY(SimpleModel::particle_var::t_div) = division_time_d.draw(gen);
    random_pool.free_state(gen);

    Uptake<UptakeDefault<typename Self::FloatType>,SimpleModel>::division(random_pool, idx, idx2, arr, buffer_arr);
  }

  KOKKOS_INLINE_FUNCTION void SimpleModel::contribution([[maybe_unused]] std::size_t idx,
                                                        std::size_t position,
                                                        double weight,
                                                        [[maybe_unused]] const SelfParticle& arr,
                                                        const MC::ContributionView& contributions)
  {
    auto access = contributions.access();
    access(position, 0) += -weight * GET_PROPERTY(SimpleModel::particle_var::phi_s); // NOLINT
  }

} // namespace Models

#endif