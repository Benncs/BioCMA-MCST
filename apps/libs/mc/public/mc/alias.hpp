#ifndef __MC__ALIAS_HPP__
#define __MC__ALIAS_HPP__

#include "Kokkos_Random.hpp"
#include <Kokkos_Core_fwd.hpp>
#include <Kokkos_ScatterView.hpp>
#include <common/traits.hpp>
#include <traits/Kokkos_IterationPatternTrait.hpp>
#include <type_traits>

// static_assert(FloatingPointType<Kokkos::Experimental::half_t>,
//               "Kokkos half_t ok");

namespace MC
{

  // NOLINTBEGIN(hicpp-avoid-c-arrays,cppcoreguidelines-avoid-c-arrays,modernize-avoid-c-arrays)
  template <uint64_t Nd, FloatingPointType F>
  using ParticlesModel = Kokkos::View<F* [Nd], Kokkos::LayoutRight>;
  template <FloatingPointType F>
  using DynParticlesModel = Kokkos::View<F**, Kokkos::LayoutRight>;
  using ComputeSpace = Kokkos::DefaultExecutionSpace;
  using HostSpace = Kokkos::DefaultHostExecutionSpace;
  // NOLINTEND(hicpp-avoid-c-arrays,cppcoreguidelines-avoid-c-arrays,modernize-avoid-c-arrays)

}; // namespace MC

using kernelMT = Kokkos::MemoryTraits<Kokkos::MemoryTraitsFlags::RandomAccess |
                                      Kokkos::MemoryTraitsFlags::Atomic |
                                      Kokkos::MemoryTraitsFlags::Restrict>;

using kernelContribution =
    Kokkos::View<float**, Kokkos::LayoutLeft, MC::ComputeSpace, kernelMT>;

namespace MC
{
  template <typename ExecSpace>
  using gen_pool_type = Kokkos::Random_XorShift1024_Pool<ExecSpace>;

  using pool_type = gen_pool_type<Kokkos::DefaultExecutionSpace>;

  using generator_type = pool_type::generator_type;

  using restrict_mt = Kokkos::MemoryTraits<Kokkos::MemoryTraitsFlags::Restrict>;
  namespace
  {

    /**
     * @brief Store particle ages
     * ages(i,0) -> hydraulic time (set to 0 when leave reactor)
     * ages(i,1) -> time since division (set to 0 when divide)
     **/
    template <typename Exec>
    using ParticleAgesBase =
        Kokkos::View<double* [2], Kokkos::LayoutLeft, Exec>; // NOLINT
  } // namespace

  enum class Status : char
  {
    Idle = 0,
    Division,
    Exit,
    Dead,
  };

  struct ContribIndexBounds
  {
    int begin;
    int end;
  };

  struct LeavingFlow;

  using ParticlePositions = Kokkos::View<uint64_t*, ComputeSpace>;
  using ParticleStatus = Kokkos::View<Status*, ComputeSpace>;
  using ParticleWeigths = Kokkos::View<double*, ComputeSpace>;
  using ParticleAges = ParticleAgesBase<ComputeSpace>;

  using ParticleSamples = Kokkos::
      View<Kokkos::Experimental::half_t**, Kokkos::LayoutRight, ComputeSpace>;

  template <typename MemorySpace>
  using ParticlePropertyViewType =
      Kokkos::View<double**, Kokkos::LayoutRight, MemorySpace>;

  using PropertySubViewtype =
      Kokkos::Subview<ParticlePropertyViewType<ComputeSpace>,
                      decltype(Kokkos::ALL),
                      std::size_t>;

  using ContributionView =
      decltype(Kokkos::Experimental::create_scatter_view(kernelContribution()));

  using KernelConcentrationType =
      Kokkos::View<const double**,
                   Kokkos::LayoutLeft,
                   ComputeSpace,
                   Kokkos::MemoryTraits<Kokkos::RandomAccess>>;
  using LocalConcentration =
      Kokkos::Subview<KernelConcentrationType, int, decltype(Kokkos::ALL)>;

  //   template <class Scalar, class ExecSpace>
  //   using VolumeView = std::enable_if_t<
  //       std::is_same_v<std::remove_const_t<Scalar>, double>,
  //       Kokkos::
  //           View<Scalar*, ExecSpace,
  //           Kokkos::MemoryTraits<Kokkos::RandomAccess>>>;

  template <class ExecSpace, bool is_const>
  using VolumeView = std::conditional_t<
      is_const,
      Kokkos::View<const double*,
                   ExecSpace,
                   Kokkos::MemoryTraits<Kokkos::RandomAccess>>,
      Kokkos::View<double*, ExecSpace>>;

  template <class ExecSpace, bool is_const>
  using DiagonalView = std::conditional_t<
      is_const,
      Kokkos::View<const double*,
                   Kokkos::LayoutLeft,
                   ExecSpace,
                   Kokkos::MemoryTraits<Kokkos::RandomAccess>>,
      Kokkos::View<double*, Kokkos::LayoutLeft, ExecSpace>>;

  template <bool is_const>
  using LeavingFlowView = std::conditional_t<
      is_const,
      Kokkos::View<const LeavingFlow*, Kokkos::SharedHostPinnedSpace>,
      Kokkos::View<LeavingFlow*, Kokkos::SharedHostPinnedSpace>>;

  template <class ExecSpace, bool is_const>
  using CumulativeProbabilityView = std::conditional_t<
      is_const,
      Kokkos::View<const double**,
                   Kokkos::LayoutRight,
                   ExecSpace,
                   Kokkos::MemoryTraits<Kokkos::RandomAccess>>,
      Kokkos::View<double**, Kokkos::LayoutRight, ExecSpace>>;

  template <class ExecSpace, bool is_const>
  using NeighborsView = std::conditional_t<
      is_const,
      Kokkos::View<const std::size_t**,
                   Kokkos::LayoutRight,
                   ExecSpace,
                   Kokkos::MemoryTraits<Kokkos::RandomAccess>>,
      Kokkos::View<std::size_t**, Kokkos::LayoutRight, ExecSpace>>;

}; // namespace MC
// FIXME

#endif
