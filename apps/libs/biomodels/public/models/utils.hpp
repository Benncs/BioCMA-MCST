#ifndef __BIO__UTILS_HPP__
#define __BIO__UTILS_HPP__
#include <Kokkos_Core.hpp>
#include <mc/prng/prng.hpp>

namespace Models
{
  static constexpr double tau_division_proba = 1e-7;

  KOKKOS_INLINE_FUNCTION bool
  check_probability_division(double d_t, double gamma, MC::KPRNG &_rng)
  {
    const double proba_div =
        (1 - Kokkos::exp(-d_t / tau_division_proba)) * gamma;

    const double x = _rng.double_uniform();
    return x < proba_div;
  }

} // namespace Models

#endif