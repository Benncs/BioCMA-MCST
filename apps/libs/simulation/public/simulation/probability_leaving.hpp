#ifndef __SIMULATION_PROBA_LEAVING_HPP__
#define __SIMULATION_PROBA_LEAVING_HPP__

#include <Kokkos_Core.hpp>
#include <common/maths.hpp>

namespace Simulation::KernelInline
{

  static constexpr bool _use_kokkos_log = true; // FIXME
  // KOKKOS_INLINE_FUNCTION bool probability_leaving(float random_number,
  //                                                 double volume,
  //                                                 double flow,
  //                                                 double dt)
  // {
  //   return (dt * flow) >
  //          (-CommonMaths::_ln<_use_kokkos_log>(random_number) * volume);
  // }

  template <typename FastSample = void>
  KOKKOS_INLINE_FUNCTION bool probability_leaving(float random_number,
                                                  double volume,
                                                  double flow,
                                                  double dt)
  {
    // Default behavior (with ln)
    return (dt * flow) >
           (-CommonMaths::_ln<_use_kokkos_log>(random_number) * volume);
  }

  // Specialization for when FastSample is provided
  template <>
  KOKKOS_INLINE_FUNCTION bool probability_leaving<void>(float random_number,
                                                        double volume,
                                                        double flow,
                                                        double dt)
  {
    // Fast version without ln
    return (dt * flow / volume) > random_number;
  }

} // namespace Simulation::KernelInline

#endif
