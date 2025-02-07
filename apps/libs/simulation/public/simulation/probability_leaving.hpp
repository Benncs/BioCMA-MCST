#ifndef __SIMULATION_PROBA_LEAVING_HPP__
#define __SIMULATION_PROBA_LEAVING_HPP__

#include <Kokkos_Core.hpp>
#include <type_traits>

static constexpr bool _use_kokkos_log = true; //FIXME

namespace Simulation::KernelInline
{

  template <bool use_kokkos_log>
  KOKKOS_INLINE_FUNCTION typename std::enable_if<!use_kokkos_log, float>::type _ln(float x)
  {
    unsigned int bx = *reinterpret_cast<unsigned int*>(&x);
    const unsigned int ex = bx >> 23;
    const signed int t = static_cast<signed int>(ex) - static_cast<signed int>(127);
    // unsigned int s = (t < 0) ? (-t) : t;
    bx = 1065353216 | (bx & 8388607);
    x = *reinterpret_cast<float*>(&bx);
    return -1.49278 + (2.11263 + (-0.729104 + 0.10969 * x) * x) * x + 0.6931471806 * t;
  }

  template <bool use_kokkos_log>
  KOKKOS_INLINE_FUNCTION typename std::enable_if<use_kokkos_log, float>::type _ln(float x)
  {
    return Kokkos::log(x);
  }

  KOKKOS_INLINE_FUNCTION bool
  probability_leaving(float random_number, double volume, double flow, double dt)
  {
    return (dt * flow) > (-_ln<_use_kokkos_log>(random_number) * volume);
  }
} // namespace Simulation::KernelInline


#endif 