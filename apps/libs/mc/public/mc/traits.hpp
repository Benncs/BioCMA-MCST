#ifndef __MODELS_TRAITS_HPP__
#define __MODELS_TRAITS_HPP__

#include <Kokkos_Core.hpp>
#include <Kokkos_Random.hpp>
#include <common/common.hpp>
#include <concepts>
#include <cstdint>


#define CHECK_MODEL(name) \
  static_assert(ModelType<name>, #name); \
  static_assert(FloatingPointType<name::FloatType>, " floatl" #name); \

namespace MC
{

  using ComputeSpace = Kokkos::DefaultExecutionSpace;
  using LocalConcentration =
      Kokkos::Subview<Kokkos::View<const double**>, int, decltype(Kokkos::ALL)>;

  using ParticlePositions = Kokkos::View<uint64_t*, ComputeSpace>;
  using pool_type = Kokkos::Random_XorShift1024_Pool<ComputeSpace>;

  enum class Status : int
  {
    Division = 0,
    Exit,
    Idle,
    Dead,
  };

  template <uint64_t Nd, FloatingPointType F> using ParticlesModel = Kokkos::View<F* [Nd]>;

}; // namespace MC

template <typename T>
concept ModelType = requires(T model,
                             T::FloatType d_t,
                             std::size_t idx,
                             std::size_t idx2,
                             const MC::ParticlesModel<T::n_var, typename T::FloatType>& arr,
                             const MC::ParticlesModel<T::n_var, typename T::FloatType>& buffer_arr,
                             const MC::LocalConcentration& c) {
  { T::n_var } -> std::convertible_to<std::size_t>;

  { T::init(idx, arr) };

  { T::update(d_t, idx, arr, c) } -> std::convertible_to<MC::Status>;

  { T::division(idx, idx2, arr, buffer_arr) };

  typename T::FloatType;
  typename T::SelfParticle;
  typename T::Self;
};

#endif