#ifndef __MODELS_TRAITS_HPP__
#define __MODELS_TRAITS_HPP__

#include <Kokkos_Core_fwd.hpp>
#include <Kokkos_Random.hpp>
#include <Kokkos_ScatterView.hpp>
#include <common/common.hpp>
#include <concepts>
#include <cstdint>

#define CHECK_MODEL(name)                                                                          \
  static_assert(ModelType<name>, #name);                                                           \
  static_assert(FloatingPointType<name::FloatType>, " floatl" #name);

namespace MC
{
  enum Status : int
  {
    Idle=0,
    Division,
    Exit,
    Dead,
  };
  
  using ComputeSpace = Kokkos::DefaultExecutionSpace;


  using ParticlePositions = Kokkos::View<uint64_t*, ComputeSpace>;
  using ParticleStatus = Kokkos::View<Status*, ComputeSpace>;
  using pool_type = Kokkos::Random_XorShift1024_Pool<ComputeSpace>;
  using ContributionView =
      Kokkos::Experimental::ScatterView<double**, Kokkos::LayoutLeft>; ///< Contribution inside the
                                                                       ///< particle's current
                                                                       ///< container
  using pool_type = Kokkos::Random_XorShift1024_Pool<ComputeSpace>;

  using KernelConcentrationType=  Kokkos::View<const double**,
                             Kokkos::LayoutLeft,
                             ComputeSpace,
                             Kokkos::MemoryTraits<Kokkos::RandomAccess>>;
  using LocalConcentration =
      Kokkos::Subview<KernelConcentrationType, int, decltype(Kokkos::ALL)>;
  

  template <uint64_t Nd, FloatingPointType F> using ParticlesModel = Kokkos::View<F* [Nd]>;

}; // namespace MC

// FIXME
using kernelContribution = Kokkos::View<double**, Kokkos::LayoutLeft, MC::ComputeSpace>;

#define GET_INDEX(size)                                                                            \
  std::size_t idx = team_handle.league_rank() * team_handle.team_size() + team_handle.team_rank(); \
  if (idx >= (size))                                                                               \
  {                                                                                                \
    return;                                                                                        \
  }

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

struct DefaultModel
{
  enum particle_var
  {
    mass = 0,
  };
  static constexpr std::size_t n_var = 1;
  static constexpr std::string_view name = "simple";
  using Self = DefaultModel;
  using FloatType = float;
  using SelfParticle = MC::ParticlesModel<Self::n_var, Self::FloatType>;

  KOKKOS_INLINE_FUNCTION static void init([[maybe_unused]] std::size_t idx,
                                          [[maybe_unused]] const SelfParticle& arr)
  {
  }

  KOKKOS_INLINE_FUNCTION static MC::Status update([[maybe_unused]] FloatType d_t,
                                                  [[maybe_unused]] std::size_t idx,
                                                  [[maybe_unused]] const SelfParticle& arr,
                                                  [[maybe_unused]] const MC::LocalConcentration& c)
  {

    return MC::Status::Idle;
  }

  KOKKOS_INLINE_FUNCTION static void division([[maybe_unused]] std::size_t idx,
                                              [[maybe_unused]] std::size_t idx2,
                                              [[maybe_unused]] const SelfParticle& arr,
                                              [[maybe_unused]] const SelfParticle& buffer_arr)
  {
  }
};

CHECK_MODEL(DefaultModel)

#endif