#ifndef __MC__ALIAS_HPP__
#define __MC__ALIAS_HPP__

#include "Kokkos_Core.hpp"
#include "decl/Kokkos_Declare_OPENMP.hpp"
#include <Kokkos_Core_fwd.hpp>
#include <Kokkos_ScatterView.hpp>
#include <common/traits.hpp>

namespace MC
{
  using ComputeSpace = Kokkos::DefaultExecutionSpace;
}

using kernelMT = Kokkos::MemoryTraits<Kokkos::MemoryTraitsFlags::RandomAccess |
                                      Kokkos::MemoryTraitsFlags::Atomic |
                                      Kokkos::MemoryTraitsFlags::Restrict>;

using kernelContribution =
    Kokkos::View<double**, Kokkos::LayoutRight, MC::ComputeSpace, kernelMT>;
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

  // template <typename Exec>
  // using GContribIndexBounds = Kokkos::View<SContribIndexBounds, Exec>; //
  // NOLINT

  // using ContribIndexBounds = GContribIndexBounds<ComputeSpace>; // NOLINT
  using restrict_mt = Kokkos::MemoryTraits<Kokkos::MemoryTraitsFlags::Restrict>;
  using ParticlePositions = Kokkos::View<uint64_t*, ComputeSpace, restrict_mt>;
  using ParticleStatus = Kokkos::View<Status*, ComputeSpace, restrict_mt>;
  using ParticleWeigths = Kokkos::View<double*, ComputeSpace, restrict_mt>;
  using ParticleAges =
      Kokkos::View<double* [2], Kokkos::LayoutLeft, ComputeSpace, restrict_mt>;

  template <typename MemorySpace>
  using ParticlePropertyViewType =
      Kokkos::View<double**, Kokkos::LayoutRight, MemorySpace>;

  using PropertySubViewtype =
      Kokkos::Subview<ParticlePropertyViewType<ComputeSpace>,
                      decltype(Kokkos::ALL),
                      std::size_t>;

  // Kernel alias
  using ContributionView =
      decltype(Kokkos::Experimental::create_scatter_view(kernelContribution()));
  ///< Contribution inside the
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

#endif
