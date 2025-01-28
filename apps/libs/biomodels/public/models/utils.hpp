#ifndef __BIO__UTILS_HPP__
#define __BIO__UTILS_HPP__
#include "mc/particles/data_holder.hpp"
#include <Kokkos_Core.hpp>
#include <mc/prng/prng.hpp>
#include <Kokkos_Random.hpp>
namespace Models
{
  static constexpr double tau_division_proba = 1e-7;

  KOKKOS_INLINE_FUNCTION bool
  check_probability_division(double /*d_t*/, double gamma, Kokkos::Random_XorShift64_Pool<Kokkos::DefaultExecutionSpace> random_pool)
  {
    // const double proba_div =
    //     (1 - Kokkos::exp(-d_t / tau_division_proba)) * gamma;
    double proba_div = gamma;
    // const double x = _rng.double_uniform();


    auto generator = random_pool.get_state();
    const double x = generator.drand(0., 1.);
    random_pool.free_state(generator);

    return x < proba_div;
  }

  
  KOKKOS_INLINE_FUNCTION void update_division_status(MC::CellStatus &status,
                                                     double d_t,
                                                     double gamma,
                                                     Kokkos::Random_XorShift64_Pool<Kokkos::DefaultExecutionSpace> rng)
  {
    status = Models::check_probability_division(d_t, gamma, rng)
                 ? MC::CellStatus::CYTOKINESIS
                 : MC::CellStatus::IDLE;


  }

  namespace GammaDivision
  {
    KOKKOS_INLINE_FUNCTION double
    threshold_linear(double lenght, double threshold, double upper_bound)
    {
      return (lenght <= threshold)
                 ? 0
                 : (lenght - threshold) / (upper_bound - threshold);
    }

  } // namespace GammaDivision


  KOKKOS_INLINE_FUNCTION bool almost_equal(double val, double val2, double tolerance=1e-8)
  {
    return Kokkos::abs(val - val2) < tolerance;
  }

  template <typename... Args> KOKKOS_INLINE_FUNCTION double min_var(Args... args)
  {
    return (Kokkos::min)({args...});
  }
  //Explicit definition to handle case where no argumet is given
  template <> KOKKOS_INLINE_FUNCTION double min_var()
  {
    return 0.;
  }

} // namespace Models

#endif
