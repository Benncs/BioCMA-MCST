#ifndef __SIMULATION__MOVE_INFO_HPP__
#define __SIMULATION__MOVE_INFO_HPP__

#include <common/common.hpp>
#include <Kokkos_Core_fwd.hpp>
#include <mc/alias.hpp>
#include <simulation/alias.hpp>

namespace Simulation::KernelInline
{

  struct LeavingFlow
  {
    std::size_t index;
    double flow;
  };

  template <typename ExecSpace> struct MoveInfo
  {
    ConstNeighborsView<ExecSpace> neighbors;
    DiagonalView<ExecSpace> diag_transition;
    CumulativeProbabilityView<ExecSpace> cumulative_probability;
    // LeavingFlowType leaving_flow;
    // LeavingFlowIndexType index_leaving_flow;
    Kokkos::View<LeavingFlow*,
                 Kokkos::SharedHostPinnedSpace,
                 Kokkos::MemoryTraits<Kokkos::MemoryTraitsFlags::Aligned |
                                      Kokkos::MemoryTraitsFlags::Restrict>>
        leaving_flow;
    Kokkos::View<double*, ExecSpace> liquid_volume;

    MoveInfo() : MoveInfo(0, 0)
    {
    }

    MoveInfo(const std::size_t n_compartments, const std::size_t n_flows)
        : diag_transition("diag_transition", n_compartments),
          leaving_flow("leaving_flow", n_flows),
          liquid_volume("liquid_volume", n_compartments)

    {
    }

    void update(std::span<const double> newliquid_volume,
                std::span<const double> diag,
                std::span<const double> proba)
    {
      KOKKOS_ASSERT(newliquid_volume.size() == this->liquid_volume.extent(0));
      const auto n_c = diag.size();
      const auto n_neighbors = proba.size() / n_c;
      KOKKOS_ASSERT(proba.size() % n_c == 0);

      Kokkos::View<const double*, HostSpace> hostli(newliquid_volume.data(),
                                                    newliquid_volume.size());
      Kokkos::deep_copy(this->liquid_volume, hostli);

      Kokkos::View<const double*,
                   Kokkos::LayoutLeft,
                   HostSpace,
                   Kokkos::MemoryTraits<Kokkos::RandomAccess>>
          _diag_transition(diag.data(), n_c);

      Kokkos::deep_copy(diag_transition, _diag_transition);

      const auto* chunk_proba = proba.data();
      CumulativeProbabilityView<HostSpace> cumproba_host(
          chunk_proba, n_c, n_neighbors);
      cumulative_probability =
          Kokkos::create_mirror_view_and_copy(ComputeSpace(), cumproba_host);
    }

    void
    set_flow(const std::size_t i, const std::size_t i_flow, const double flow)
    {
      leaving_flow(i) = {i_flow, flow};
      // Kokkos::deep_copy(leaving_flow(i), LeavingFlow{i_flow, flow});
    }
  };
} // namespace Simulation::KernelInline
#endif
