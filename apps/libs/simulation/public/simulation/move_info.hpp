#ifndef __SIMULATION__MOVE_INFO_HPP__
#define __SIMULATION__MOVE_INFO_HPP__

#include <Kokkos_Core_fwd.hpp>
#include <cma_utils/cache_hydro_state.hpp>
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

    void
    set_flow(const std::size_t i, const std::size_t i_flow, const double flow)
    {
      leaving_flow(i) = {i_flow, flow};
      // Kokkos::deep_copy(leaving_flow(i), LeavingFlow{i_flow, flow});
    }
  };
} // namespace Simulation::KernelInline
#endif
