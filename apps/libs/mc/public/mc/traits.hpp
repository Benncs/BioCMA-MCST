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
#include <optional>
#include <type_traits>

/** @brief Utility for compile time array concatenation  */
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
  template <uint64_t Nd, FloatingPointType F>
  using ParticlesModel = Kokkos::View<F* [Nd], Kokkos::LayoutRight>;
  template <FloatingPointType F>
  using DynParticlesModel = Kokkos::View<F**, Kokkos::LayoutRight>;
  // NOLINTEND(hicpp-avoid-c-arrays,cppcoreguidelines-avoid-c-arrays,modernize-avoid-c-arrays)

}; // namespace MC

template <typename T>
concept ConfigurableInit = requires(T model,
                                    const std::size_t size,
                                    const MC::KPRNG::pool_type& random_pool,
                                    std::size_t idx,
                                    const typename T::SelfParticle& arr,
                                    const T::Config& config) {
  { model.init(random_pool, idx, arr, config) } -> std::same_as<void>;

  { model.get_config(size) } -> std::same_as<typename T::Config>;
};

template <typename T>
concept NonConfigurableInit = requires(T model,
                                       const MC::KPRNG::pool_type& random_pool,
                                       std::size_t idx,
                                       const typename T::SelfParticle& arr) {
  { model.init(random_pool, idx, arr) } -> std::same_as<void>;
};
using NonConfigType = std::nullopt_t;
/**
  @brief Concept to define a correct Model

  Minimum model definition should validate this contract, the behaviour behind
  each function can vary but the order and signature should be the same
 */
template <typename T, typename ViewType>
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
                                   const MC::KPRNG::pool_type& random_pool,
                                   const T::Config& config) {
  {
    T::n_var
  } -> std::convertible_to<std::size_t>; ///< A model should declare the number
                                         ///< of internal variable
  typename T::FloatType; ///< Type used internally by model to declare internal
                         ///< floating point values
  typename T::SelfParticle; ///< Equivalent to MC::ParticlesModel<Self::n_var,
                            ///< Self::FloatType>
  typename T::Self;         ///< Model typename
  typename T::Config;       ///< Model typename

  // Check if the model is configurable
  // requires(std::is_same_v<typename T::Config, std::nullopt_t> ?
  // NonConfigurableInit<T>
  //                                                             :
  //                                                             ConfigurableInit<T>);

  requires ConfigurableInit<T> ||
               (std::is_same_v<typename T::Config, NonConfigType> &&
                NonConfigurableInit<T>);

  // {
  //   T::init(random_pool, idx, arr, config)
  // } -> std::same_as<void>; ///< Main init function applied to each MC
  // particle at the begin of the
  //                          ///< simulation

  {
    T::mass(idx, arr)
  } -> std::same_as<double>; ///< Return the individual mass of particle

  {
    T::update(random_pool, d_t, idx, arr, c)
  } -> std::convertible_to<MC::Status>; ///< Update state of MC particle

  // {
  //   T::contribution(idx, position, weight, arr, contributions)
  // } -> std::same_as<void>; ///< Get the individual contribution for the MC
  // particle

  { T::get_bounds() } -> std::same_as<MC::ContribIndexBounds>;

  {
    T::division(random_pool, idx, idx2, arr, buffer_arr)
  } -> std::same_as<void>; ///< Perform the internal property redistribution
                           ///< after division

  requires FloatingPointType<typename T::FloatType>;
};

/** @brief  SFNIAE way to declare a model with number of internal properties
 * known at compile time
 */
template <typename T>
concept FixedModelType =
    CommonModelType<T, MC::ParticlesModel<T::n_var, typename T::FloatType>>;

/** @brief  SFNIAE wau to declare a model with number of internal properties not
known at compile time Alows to properly define User Defined Model but should be
avoid in other cases
 */
template <typename T>
concept DynModelType =
    CommonModelType<T, MC::DynParticlesModel<typename T::FloatType>>;

/** @brief  Model type
 */
template <typename T>
concept ModelType = DynModelType<T> || FixedModelType<T>;

template <typename T>
concept NonConfigurableModel = ModelType<T> && NonConfigurableInit<T>;

template <typename T>
concept ConfigurableModel = ModelType<T> && ConfigurableInit<T>;

/** @brief SFNIAE way to check whether model allow internal value saving or not
 */
template <std::size_t n, typename T>
concept _HasExportProperties = requires(const T obj) {
  { T::names() } -> std::convertible_to<std::array<std::string_view, n>>;
};

/** @brief SFNIAE way to check whether model allow all value saving*/
template <typename T>
concept HasExportPropertiesFull = FixedModelType<T> && requires(const T obj) {
  { T::names() } -> std::convertible_to<std::array<std::string_view, T::n_var>>;
};

/** @brief SFNIAE way to check whether model allow partial value saving */
template <typename T>
concept HasExportPropertiesPartial = ModelType<T> && requires(const T obj) {
  { T::names() } -> std::convertible_to<std::vector<std::string_view>>;
  { T::get_number() } -> std::convertible_to<std::vector<std::size_t>>;
};

/** @brief Model that can export properties */
template <typename T>
concept HasExportProperties =
    HasExportPropertiesFull<T> || HasExportPropertiesPartial<T>;

// Helper to detect if `uniform_weight` exists as a type alias (using `using`
// keyword)
template <typename T, typename = void>
struct has_uniform_weight : std::false_type
{
};

template <typename T>
struct has_uniform_weight<T, std::void_t<typename T::uniform_weight>>
    : std::true_type
{
};

/** @brief Concept to check if a model type has `uniform_weight`*/
template <typename T>
concept ConstWeightModelType = ModelType<T> && has_uniform_weight<T>::value;

/** @brief Concept to check if a model type has `uniform_weight`*/
template <typename T>
concept PreInitModel = ModelType<T> && requires(T model) { T::preinit(); };

#endif
