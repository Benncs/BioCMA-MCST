#ifndef __MC__ALIAS_HPP__
#define __MC__ALIAS_HPP__

#include <Kokkos_Core_fwd.hpp>
#include <Kokkos_Random.hpp>
#include <Kokkos_ScatterView.hpp>
#include <common/traits.hpp>
#include <decl/Kokkos_Declare_OPENMP.hpp>
#include <traits/Kokkos_IterationPatternTrait.hpp>
#include <type_traits>

// static_assert(FloatingPointType<Kokkos::Experimental::half_t>,
//               "Kokkos half_t ok");

namespace MC
{
  using ComputeSpace = Kokkos::DefaultExecutionSpace;
  using HostSpace = Kokkos::DefaultHostExecutionSpace;

  using ComputeLayout = ComputeSpace::array_layout;

  // NOLINTBEGIN(hicpp-avoid-c-arrays,cppcoreguidelines-avoid-c-arrays,modernize-avoid-c-arrays)
  template <uint64_t Nd, FloatingPointType F>
  //   using ParticlesModel = Kokkos::View<
  //       F* [Nd],
  //       Kokkos::LayoutRight,
  //       ComputeSpace,
  //       Kokkos::MemoryTraits<Kokkos::MemoryTraitsFlags::Restrict> >;

  //   template <uint64_t Nc, FloatingPointType F>
  //   using ParticlesContribs = Kokkos::View<
  //       F* [Nc],
  //       Kokkos::LayoutRight,
  //       ComputeSpace,
  //       Kokkos::MemoryTraits<Kokkos::MemoryTraitsFlags::Restrict
  //                            | Kokkos::MemoryTraitsFlags::Aligned> >;

  //   template <FloatingPointType F>
  //   using DynParticlesModel = Kokkos::View<
  //       F**,
  //       Kokkos::LayoutRight,
  //       ComputeSpace,
  //       Kokkos::MemoryTraits<Kokkos::MemoryTraitsFlags::Restrict> >;

  //   template <FloatingPointType F>
  //   using DynParticlesContribs = Kokkos::View<
  //       F**,
  //       Kokkos::LayoutRight,
  //       ComputeSpace,
  //       Kokkos::MemoryTraits<Kokkos::MemoryTraitsFlags::Restrict> >;

  using ParticlesModel
      = Kokkos::View<F* [Nd],
                     ComputeSpace::array_layout,
                     ComputeSpace,
                     Kokkos::MemoryTraits<Kokkos::MemoryTraitsFlags::Restrict>>;

  template <uint64_t Nc, FloatingPointType F>
  using ParticlesContribs = Kokkos::View<
      F* [Nc],
      ComputeSpace::array_layout,
      ComputeSpace,
      Kokkos::MemoryTraits<Kokkos::MemoryTraitsFlags::Restrict
                           | Kokkos::MemoryTraitsFlags::Aligned>>;

  template <FloatingPointType F>
  using DynParticlesModel
      = Kokkos::View<F**,
                     ComputeSpace::array_layout,
                     ComputeSpace,
                     Kokkos::MemoryTraits<Kokkos::MemoryTraitsFlags::Restrict>>;

  template <FloatingPointType F>
  using DynParticlesContribs
      = Kokkos::View<F**,
                     ComputeSpace::array_layout,
                     ComputeSpace,
                     Kokkos::MemoryTraits<Kokkos::MemoryTraitsFlags::Restrict>>;

  // NOLINTEND(hicpp-avoid-c-arrays,cppcoreguidelines-avoid-c-arrays,modernize-avoid-c-arrays)

}; // namespace MC

using kernelMT = Kokkos::MemoryTraits<Kokkos::MemoryTraitsFlags::RandomAccess
                                      | Kokkos::MemoryTraitsFlags::Atomic
                                      | Kokkos::MemoryTraitsFlags::Restrict>;

using kernelContribution
    = Kokkos::View<float**, Kokkos::LayoutLeft, MC::ComputeSpace, kernelMT>;

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
    // template <typename Exec>
    // using ParticleAgesBase = Kokkos::
    //     View<double* [2], typename Exec::array_layout, Exec>; // NOLINT

    // Force layout left because iteration on age(i,0)->age(i,1) is never done
    // in compute kernel
    template <typename Exec>
    using ParticleAgesBase
        = Kokkos::View<float* [2], Kokkos::LayoutLeft, Exec>; // NOLINT
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

  // using ParticleSamples = Kokkos::
  //     View<Kokkos::Experimental::half_t**, Kokkos::LayoutRight,
  //     ComputeSpace>;

  // V1.0-rc: LayoutRight seems to lead to better perf than Left for GPU
  using ParticleSamples
      = Kokkos::View<float**, Kokkos::LayoutRight, ComputeSpace>;

  template <typename MemorySpace>
  using ParticlePropertyViewType
      = Kokkos::View<double**, Kokkos::LayoutRight, MemorySpace>;

  using PropertySubViewtype
      = Kokkos::Subview<ParticlePropertyViewType<ComputeSpace>,
                        decltype(Kokkos::ALL),
                        std::size_t>;

  using ContributionView = decltype(Kokkos::Experimental::create_scatter_view(
      kernelContribution()));

  // //
  // using ContributionView = class Kokkos::Experimental::ScatterView<
  //     float**,
  //     Kokkos::LayoutLeft,
  //     ComputeSpace,
  //     Kokkos::Experimental::ScatterSum,
  //     Kokkos::Experimental::ScatterNonDuplicated,
  //     Kokkos::Experimental::ScatterAtomic>;

  using KernelConcentrationType
      = Kokkos::View<const double**,
                     Kokkos::LayoutLeft,
                     ComputeSpace,
                     Kokkos::MemoryTraits<Kokkos::RandomAccess>>;

  using LocalConcentration = KernelConcentrationType;
  //   using LocalConcentration
  //       = Kokkos::Subview<KernelConcentrationType, int,
  //       decltype(Kokkos::ALL)>;

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

  // TODO Use execspace layout
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
