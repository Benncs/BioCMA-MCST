#ifndef __MODELS_TRAITS_HPP__
#define __MODELS_TRAITS_HPP__

#include <Kokkos_Core_fwd.hpp>
#include <Kokkos_Random.hpp>
#include <Kokkos_ScatterView.hpp>
#include <common/common.hpp>
#include <common/traits.hpp>
#include <concepts>
#include <cstdint>
#include <mc/alias.hpp>
#include <mc/macros.hpp>
#include <mc/prng/prng.hpp>
#include <type_traits>



template <std::size_t N1, std::size_t N2>
constexpr std::array<std::string_view, N1 + N2>
concat_arrays(const std::array<std::string_view, N1>& arr1,
              const std::array<std::string_view, N2>& arr2)
{
  std::array<std::string_view, N1 + N2> result;
  std::copy(arr1.begin(), arr1.end(), result.begin());
  std::copy(arr2.begin(), arr2.end(), result.begin() + N1);
  return result;
}

namespace MC
{

  // NOLINTBEGIN(hicpp-avoid-c-arrays,cppcoreguidelines-avoid-c-arrays,modernize-avoid-c-arrays)
  template <uint64_t Nd, FloatingPointType F> using ParticlesModel = Kokkos::View<F* [Nd]>;
  template <FloatingPointType F> using DynParticlesModel = Kokkos::View<F**>;
  // NOLINTEND(hicpp-avoid-c-arrays,cppcoreguidelines-avoid-c-arrays,modernize-avoid-c-arrays)

}; // namespace MC

/**
  @brief Concept to define a correct Model
 */
template <typename T,typename ViewType>
concept CommonModelType = requires(T model,
                                T::FloatType d_t,
                                std::size_t idx,
                                std::size_t idx2,
                                double weight,
                                const T::SelfParticle& arr,
                                const T::SelfParticle& buffer_arr,
                                const MC::LocalConcentration& c,
                                std::size_t position,
                                const MC::ContributionView& contributions,
                                const MC::KPRNG::pool_type& random_pool) {
  { T::n_var } -> std::convertible_to<std::size_t>;
  typename T::FloatType;
  typename T::SelfParticle;
  typename T::Self;
  requires FloatingPointType<typename T::FloatType>;
  { T::n_var } -> std::convertible_to<std::size_t>;

  { T::init(random_pool, idx, arr) } -> std::same_as<void>;

  { T::mass(idx, arr) } -> std::same_as<double>;

  { T::update(random_pool, d_t, idx, arr, c) } -> std::convertible_to<MC::Status>;

  { T::contribution(idx, position, weight, arr, contributions) } -> std::same_as<void>;

  { T::division(random_pool, idx, idx2, arr, buffer_arr) } -> std::same_as<void>;
};


template <typename T>
concept DynModelType = CommonModelType<T,MC::DynParticlesModel<typename T::FloatType>>;

template <typename T>
concept FixedModelType = CommonModelType<T,MC::ParticlesModel<T::n_var, typename T::FloatType>>;


template <typename T>
concept ModelType = DynModelType<T> ||  FixedModelType<T>;


template <std::size_t n, typename T>
concept _HasExportProperties = requires(const T obj) {
  { T::names() } -> std::convertible_to<std::array<std::string_view, n>>;
};

template <typename T>
concept HasExportPropertiesFull = FixedModelType<T> && requires(const T obj) {
  { T::names() } -> std::convertible_to<std::array<std::string_view, T::n_var>>;
};

template <typename T>
concept HasExportPropertiesPartial = ModelType<T> && requires(const T obj) {
  { T::names() } -> std::convertible_to<std::vector<std::string_view>>;
  { T::get_number() } -> std::convertible_to<std::vector<std::size_t>>;
};

template <typename T>
concept HasExportProperties = HasExportPropertiesFull<T> || HasExportPropertiesPartial<T>;

// Helper to detect if `uniform_weight` exists as a type alias (using `using` keyword)
template <typename T, typename = void> struct has_uniform_weight : std::false_type
{
};

template <typename T>
struct has_uniform_weight<T, std::void_t<typename T::uniform_weight>> : std::true_type
{
};

// Concept to check if a model type has `uniform_weight`
template <typename T>
concept ConstWeightModelType = ModelType<T> && has_uniform_weight<T>::value;

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

CHECK_MODEL(DynamicDefaultModel)

// namespace MC
// {
//   // TODO MOVE

//   struct TagDetector
//   {
//     KOKKOS_FUNCTION void
//     operator()(const Kokkos::TeamPolicy<ComputeSpace>::member_type& team_handle,
//                int& dead_count) const
//     {
//       (void)team_handle;
//       (void)dead_count;
//     }
//     TagDetector() = default;
//   };
//   Kokkos::TeamPolicy<ComputeSpace> inline get_policy_auto(std::size_t range)
//   {

//     Kokkos::TeamPolicy<ComputeSpace> _policy;

//     int recommended_team_size =
//         _policy.team_size_recommended(TagDetector(), Kokkos::ParallelForTag());
//     int league_size = (static_cast<int>(range) + recommended_team_size - 1) /
//     recommended_team_size;

//     _policy = Kokkos::TeamPolicy<ComputeSpace>(league_size, recommended_team_size);

//     return _policy;
//   }

//   template <typename FunctorType>
//   Kokkos::TeamPolicy<ComputeSpace>
//   get_policy(FunctorType& f, std::size_t range, bool reduce = false)
//   {

//     // Kokkos::TeamPolicy<ComputeSpace> _policy;
//     // int recommended_team_size = (reduce)
//     //                                 ? _policy.team_size_recommended(f,
//     //                                 Kokkos::ParallelReduceTag()) :
//     //                                 _policy.team_size_recommended(f,
//     Kokkos::ParallelForTag());

//     // int league_size = (static_cast<int>(range) + recommended_team_size - 1) /
//     // recommended_team_size;

//     // return {league_size, recommended_team_size};

//     Kokkos::TeamPolicy<ComputeSpace> _policy;

//     int recommended_team_size =
//         _policy.team_size_recommended(TagDetector(), Kokkos::ParallelForTag());
//     int league_size = (static_cast<int>(range) + recommended_team_size - 1) /
//     recommended_team_size;

//     _policy = Kokkos::TeamPolicy<ComputeSpace>(league_size, recommended_team_size);

//     return _policy;
//   }
// } // namespace MC

#endif