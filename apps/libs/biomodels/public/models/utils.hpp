#ifndef __BIO__UTILS_HPP__
#define __BIO__UTILS_HPP__

#include "Kokkos_MathematicalConstants.hpp"
#include "common/traits.hpp"
#include <mc/particles/data_holder.hpp>
#include <Kokkos_Core.hpp>
#include <Kokkos_Random.hpp>
#include <mc/prng/prng.hpp>

/**
  @brief Models definition
 */
namespace Models
{

  namespace MolarMass
  {
    constexpr double glucose = 180e-3;
    constexpr double dioxygen = 32e-3;
    constexpr double acetate = 59e-3;
    constexpr double co2 = 44e-3;
    constexpr double X = 113.1e-3;
  }; // namespace MolarMass


  template<FloatingPointType F>
  KOKKOS_INLINE_FUNCTION consteval F c_linear_density(F rho,F d)
  {
    return rho * F(Kokkos::numbers::pi) * d * d / F(4.);
  }

  static constexpr double tau_division_proba = 1e-7;

  KOKKOS_INLINE_FUNCTION bool check_probability_division(
      double d_t,
      double gamma,
      Kokkos::Random_XorShift64_Pool<Kokkos::DefaultExecutionSpace> random_pool)
  {
    (void)d_t;
    // const double proba_div =
    //     (1 - Kokkos::exp(-d_t / tau_division_proba)) * gamma;
 
    auto generator = random_pool.get_state();
    const double x = generator.drand(0., 1.);
    random_pool.free_state(generator);

    return x < gamma;
  }

  KOKKOS_INLINE_FUNCTION void
  update_division_status(MC::CellStatus& status,
                         double d_t,
                         double gamma,
                         Kokkos::Random_XorShift64_Pool<Kokkos::DefaultExecutionSpace> rng)
  {
    status = Models::check_probability_division(d_t, gamma, rng) ? MC::CellStatus::CYTOKINESIS
                                                                 : MC::CellStatus::IDLE;
  }

  namespace GammaDivision
  {
    KOKKOS_INLINE_FUNCTION double
    threshold_linear(double lenght, double threshold, double upper_bound)
    {
      return (lenght <= threshold) ? 0 : (lenght - threshold) / (upper_bound - threshold);
    }

    KOKKOS_INLINE_FUNCTION constexpr float logistic(float x, float xmax, float alpha)
    {
      // auto blna = alpha * std::log(xmax);
      // return 1 / (1 + std::exp(-alpha * std::log(x) + blna));

      // return x<xmax?1. / (1 + std::exp(alpha * std::log((x-xmax)/xmax ))):1.;

      const auto z = x<xmax? std::pow(x/(xmax-x),alpha):1.;
      return z/(1.+z);
    }

  } // namespace GammaDivision

  KOKKOS_INLINE_FUNCTION bool almost_equal(double val, double val2, double tolerance = 1e-8)
  {
    return Kokkos::abs(val - val2) < tolerance;
  }

  template <typename... Args> KOKKOS_INLINE_FUNCTION double min_var(Args... args)
  {
    return (Kokkos::min)({args...});
  }
  template <> KOKKOS_INLINE_FUNCTION double min_var()
  {
    return 0.;
  }

} // namespace Models

#endif
