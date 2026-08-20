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

  // template <typename T>
  // KOKKOS_INLINE_FUNCTION bool
  // bernouilli_proba(T volume, T flow, T dt)
  // {
  //   const auto _lambda = dt * flow / volume;
  //   return T{ 1 } - Kokkos::exp(-_lambda);
  // }

  template <typename T, typename FastSample = precision_tag>
  KOKKOS_INLINE_FUNCTION bool
  probability_leaving(const T random_number, const double lambda)
  {
    return lambda > -CommonMaths::_ln<_use_kokkos_log>(random_number);
  }

  template <typename T, typename FastSample = precision_tag>
  KOKKOS_INLINE_FUNCTION bool
  probability_leaving(T random_number, double volume, double flow, double dt)
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
  probability_leaving<float, fast_tag>(float random_number,
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

  // Specialization for when FastSample is provided
  template <>
  KOKKOS_INLINE_FUNCTION bool
  probability_leaving<double, fast_tag>(double random_number,
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
