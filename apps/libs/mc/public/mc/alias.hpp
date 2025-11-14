#ifndef __MC__ALIAS_HPP__
#define __MC__ALIAS_HPP__

#include <Kokkos_Core_fwd.hpp>
#include <Kokkos_ScatterView.hpp>
#include <common/traits.hpp>

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

  using ParticlePositions = Kokkos::View<uint64_t*, ComputeSpace>;
  using ParticleStatus = Kokkos::View<Status*, ComputeSpace>;
  using ParticleWeigths = Kokkos::View<double*, ComputeSpace>;
  using ParticleAges = ParticleAgesBase<ComputeSpace>;

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

}; // namespace MC
// FIXME
template <typename Space>
using ConstNeighborsView =
    Kokkos::View<const std::size_t**,
                 Kokkos::LayoutRight,
                 Space,
                 Kokkos::MemoryTraits<Kokkos::RandomAccess>>;

#endif
