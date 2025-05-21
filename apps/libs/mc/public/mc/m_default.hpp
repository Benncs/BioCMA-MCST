#ifndef __DEFAULT_MODEL_IMPLEMENTATION_HPP__
#define __DEFAULT_MODEL_IMPLEMENTATION_HPP__

#include <mc/traits.hpp>
#include <optional>


/***
DEFAULT MODEL IMPLEMENTATION
***/

struct DefaultModel
{
  enum class particle_var : int
  {
    mass = 0,
  };
  static constexpr std::size_t n_var = 1;
  
  static constexpr std::string_view name = "simple";
  using uniform_weight = std::true_type; // Using type alias
  using Self = DefaultModel;
  using FloatType = float;
  using SelfParticle = MC::ParticlesModel<Self::n_var, Self::FloatType>;
  using Config = std::nullopt_t;
  static constexpr bool uniform_weigth = false;

  KOKKOS_INLINE_FUNCTION static void init([[maybe_unused]] const MC::KPRNG::pool_type& random_pool,
                                          [[maybe_unused]] std::size_t idx,
                                          [[maybe_unused]] const SelfParticle& arr)
  {
  }

  KOKKOS_INLINE_FUNCTION static double mass([[maybe_unused]] std::size_t idx,
                                            [[maybe_unused]] const SelfParticle& arr)
  {
    return 1.;
  }

  KOKKOS_INLINE_FUNCTION static MC::Status
  update([[maybe_unused]] const MC::KPRNG::pool_type& random_pool,
         [[maybe_unused]] FloatType d_t,
         [[maybe_unused]] std::size_t idx,
         [[maybe_unused]] const SelfParticle& arr,
         [[maybe_unused]] const MC::LocalConcentration& c)
  {

    return MC::Status::Idle;
  }

  KOKKOS_INLINE_FUNCTION static void
  division([[maybe_unused]] const MC::KPRNG::pool_type& random_pool,
           [[maybe_unused]] std::size_t idx,
           [[maybe_unused]] std::size_t idx2,
           [[maybe_unused]] const SelfParticle& arr,
           [[maybe_unused]] const SelfParticle& buffer_arr)
  {
  }

  KOKKOS_INLINE_FUNCTION static void
  contribution([[maybe_unused]] std::size_t idx,
               [[maybe_unused]] std::size_t position,
               [[maybe_unused]] double weight,
               [[maybe_unused]] const SelfParticle& arr,
               [[maybe_unused]] const MC::ContributionView& contributions)
  {
  }
};

CHECK_MODEL(DefaultModel)

struct DynamicDefaultModel
{
  enum class particle_var : int
  {
    mass = 0,
  };
  static constexpr std::size_t n_var = 1;
  static constexpr std::string_view name = "simple";
  using uniform_weight = std::true_type; // Using type alias
  using Self = DynamicDefaultModel;
  using FloatType = float;
  using SelfParticle = MC::DynParticlesModel<FloatType>;
  using Config = float;
  static constexpr bool uniform_weigth = false;

  KOKKOS_INLINE_FUNCTION static void init([[maybe_unused]] const MC::KPRNG::pool_type& random_pool,
                                          [[maybe_unused]] std::size_t idx,
                                          [[maybe_unused]] const SelfParticle& arr,[[maybe_unused]] const Config& config)
  {
  }

  KOKKOS_INLINE_FUNCTION static double mass([[maybe_unused]] std::size_t idx,
                                            [[maybe_unused]] const SelfParticle& arr)
  {
    return 1.;
  }

  KOKKOS_INLINE_FUNCTION static MC::Status
  update([[maybe_unused]] const MC::KPRNG::pool_type& random_pool,
         [[maybe_unused]] FloatType d_t,
         [[maybe_unused]] std::size_t idx,
         [[maybe_unused]] const SelfParticle& arr,
         [[maybe_unused]] const MC::LocalConcentration& c)
  {

    return MC::Status::Idle;
  }

  KOKKOS_INLINE_FUNCTION static void
  division([[maybe_unused]] const MC::KPRNG::pool_type& random_pool,
           [[maybe_unused]] std::size_t idx,
           [[maybe_unused]] std::size_t idx2,
           [[maybe_unused]] const SelfParticle& arr,
           [[maybe_unused]] const SelfParticle& buffer_arr)
  {
  }

  KOKKOS_INLINE_FUNCTION static void
  contribution([[maybe_unused]] std::size_t idx,
               [[maybe_unused]] std::size_t position,
               [[maybe_unused]] double weight,
               [[maybe_unused]] const SelfParticle& arr,
               [[maybe_unused]] const MC::ContributionView& contributions)
  {
  }
};

CHECK_MODEL(DynamicDefaultModel)



#endif 