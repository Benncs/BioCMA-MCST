#ifndef __BIOMC_COMMON_MATHS__
#define __BIOMC_COMMON_MATHS__

#include <Kokkos_Core.hpp>

namespace CommonMaths
{

  template <bool use_kokkos_log>
  KOKKOS_INLINE_FUNCTION float _ln(float x)
    requires(!use_kokkos_log)
  {
    // NOLINTBEGIN
    unsigned int bx = *reinterpret_cast<unsigned int*>(&x);
    const unsigned int ex = bx >> 23;
    const signed int t =
        static_cast<signed int>(ex) - static_cast<signed int>(127);
    // unsigned int s = (t < 0) ? (-t) : t;
    bx = 1065353216 | (bx & 8388607);
    x = *reinterpret_cast<float*>(&bx);
    return -1.49278 + (2.11263 + (-0.729104 + 0.10969 * x) * x) * x +
           0.6931471806 * t;
    // NOLINTEND
  }
  template <bool use_kokkos_log>
  KOKKOS_INLINE_FUNCTION float _ln(float x)
    requires(use_kokkos_log)
  {
    return Kokkos::log(x);
  }

} // namespace CommonMaths
#endif
