#ifndef __BIO__UTILS_HPP__
#define __BIO__UTILS_HPP__
#include "mc/particles/data_holder.hpp"
#include <Kokkos_Core.hpp>
#include <mc/prng/prng.hpp>

namespace Models
{
  static constexpr double tau_division_proba = 1e-7;

  KOKKOS_INLINE_FUNCTION bool
  check_probability_division(double /*d_t*/, double gamma, MC::KPRNG &_rng)
  {
    // const double proba_div =
    //     (1 - Kokkos::exp(-d_t / tau_division_proba)) * gamma;
    double proba_div = gamma;
    const double x = _rng.double_uniform();

    return x < proba_div;
  }

  template <typename F>
  KOKKOS_INLINE_FUNCTION void update_division_status(MC::CellStatus &status,
                                                     double d_t,
                                                     F predicate,
                                                     MC::KPRNG &rng)
  {
    status = Models::check_probability_division(d_t, predicate, rng)
                 ? MC::CellStatus::CYTOKINESIS
                 : MC::CellStatus::IDLE;

  }

  namespace GammaDivision
  {
    KOKKOS_INLINE_FUNCTION double
    threshold_linear(double lenght, double threshold, double upper_bound)
    {
      return (lenght < threshold)
                 ? 0
                 : (lenght - threshold) / (upper_bound - threshold);
    }

  } // namespace GammaDivision

} // namespace Models

#endif