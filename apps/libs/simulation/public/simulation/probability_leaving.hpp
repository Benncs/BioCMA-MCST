#ifndef __SIMULATION_PROBA_LEAVING_HPP__
#define __SIMULATION_PROBA_LEAVING_HPP__

#include <Kokkos_Assert.hpp>
#include <Kokkos_Core.hpp>
#include <common/maths.hpp>

namespace Simulation::KernelInline
{

  using fast_tag = void;
  using precision_tag = int;

  static constexpr bool _use_kokkos_log = true; // FIXME

  template <typename FastSample = precision_tag>
  KOKKOS_INLINE_FUNCTION bool
  probability_leaving(float random_number,
                      double volume,
                      double flow,
                      double dt)
  {
    KOKKOS_ASSERT(random_number >= 0. && random_number <= 1.);
    KOKKOS_ASSERT(volume >= 0.);
    KOKKOS_ASSERT(flow >= 0.);
    KOKKOS_ASSERT(dt >= 0.);
    // Default behavior (with ln)
    return (dt * flow)
           > (-CommonMaths::_ln<_use_kokkos_log>(random_number) * volume);
  }

  // Specialization for when FastSample is provided
  template <>
  KOKKOS_INLINE_FUNCTION bool
  probability_leaving<fast_tag>(float random_number,
                                double volume,
                                double flow,
                                double dt)
  {
    KOKKOS_ASSERT(random_number >= 0. && random_number <= 1.);
    KOKKOS_ASSERT(volume >= 0.);
    KOKKOS_ASSERT(flow >= 0.);
    KOKKOS_ASSERT(dt >= 0.);
    // Fast version without ln
    return (dt * flow / volume) > random_number;
  }

} // namespace Simulation::KernelInline

#endif
