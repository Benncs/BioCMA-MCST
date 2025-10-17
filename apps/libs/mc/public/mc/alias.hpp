#ifndef __MC__ALIAS_HPP__
#define __MC__ALIAS_HPP__

#include <Kokkos_Core_fwd.hpp>
#include <Kokkos_ScatterView.hpp>
#include <common/traits.hpp>

namespace MC
{
  enum class Status : int
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

  using ComputeSpace = Kokkos::DefaultExecutionSpace;

  // template <typename Exec>
  // using GContribIndexBounds = Kokkos::View<SContribIndexBounds, Exec>; //
  // NOLINT

  // using ContribIndexBounds = GContribIndexBounds<ComputeSpace>; // NOLINT

  using ParticlePositions = Kokkos::View<uint64_t*, ComputeSpace>;
  using ParticleStatus = Kokkos::View<Status*, ComputeSpace>;
  using ParticleWeigths = Kokkos::View<double*, ComputeSpace>;
  using ParticleAges =
      Kokkos::View<double* [2], Kokkos::LayoutLeft, ComputeSpace>;

  template <typename MemorySpace>
  using ParticlePropertyViewType =
      Kokkos::View<double**, Kokkos::LayoutRight, MemorySpace>;

  using PropertySubViewtype =
      Kokkos::Subview<ParticlePropertyViewType<ComputeSpace>,
                      decltype(Kokkos::ALL),
                      std::size_t>;

  // Kernel alias
  using ContributionView = Kokkos::Experimental::
      ScatterView<double**, Kokkos::LayoutRight>; ///< Contribution inside the
                                                  ///< particle's current
                                                  ///< container
  // using access_type = decltype(std::declval<ContributionView>().access());

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
using kernelContribution =
    Kokkos::View<double**, Kokkos::LayoutRight, MC::ComputeSpace>;

#endif
