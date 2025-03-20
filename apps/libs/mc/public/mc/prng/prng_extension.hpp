#ifndef __PRNG_EXTENSION_HPP__
#define __PRNG_EXTENSION_HPP__

#include <Kokkos_Core.hpp>
#include <Kokkos_MathematicalConstants.hpp>
#include <Kokkos_Random.hpp>
#include <cmath>
#include <common/traits.hpp>

namespace MC::Distributions
{
  /**
   @brief Concept for probability distribution laws.

   The `ProbabilityLaw` concept defines the requirements for a type `T`
   to model a probability distribution with floating-point computations.
   It ensures that the type provides methods for drawing random samples
   and computing common statistical measures such as mean, variance, and skewness.

   @tparam T The type representing a probability distribution.
   @tparam F The floating-point type used for computations.
   @tparam DeviceType The execution device type for random number generation.

   Requirements:
   - `T` must support sampling using a Kokkos random number generator.
   - `T` must provide methods to compute statistical properties:
       - `mean()`: Returns the expected value (mean) of the distribution.
       - `var()`: Returns the variance of the distribution.
       - `skewness()`: Returns the skewness, measuring asymmetry.

   Example usage:
   @code
   struct NormalDistribution {
       double mean() const { return mu; }
       double var() const { return sigma * sigma; }
       double skewness() const { return 0.0; }
       double draw(Kokkos::Random_XorShift1024<DeviceType>& gen) {return 0.;    }
     };
     static_assert(ProbabilityLaw<NormalDistribution, double, DeviceType>);
     @endcode
   */
  template <typename T, typename F, class DeviceType>
  concept ProbabilityLaw =
      FloatingPointType<F> && requires(const T& obj, Kokkos::Random_XorShift1024<DeviceType>& gen) {
        { obj.draw(gen) } -> std::same_as<F>;
        { obj.mean() } -> std::same_as<F>;
        { obj.var() } -> std::same_as<F>;
        { obj.skewness() } -> std::same_as<F>;
      };

  /**
  @brief Computes an approximation of the inverse error function.

  This function approximates the inverse error function `erfinv(x)`
  using Winitzki’s method (from A. Soranzo, E. Epure), which provides a simple and computationally
  efficient formula based on logarithms and square roots.

  @tparam F Floating-point type (must satisfy `FloatingPointType`).
  @param x Input value in the range [-1, 1].
  @return Approximate value of `erfinv(x)`.`.

  @note This implementation is not highly optimized for GPUs.
        For a more accurate and efficient implementation,
        see: https://people.maths.ox.ac.uk/gilesm/codes/erfinv/gems.pdf. See also bramowitz and
  Stegun method (Handbook of Mathematical Functions, formula 7.1.26)

  @warning This approximation is not valid for extreme values |x| close to 1.
*/
  template <FloatingPointType F> KOKKOS_INLINE_FUNCTION F erfinv(F x)
  {

    // Use the Winitzki’s method to calculate get an approached expression of erf(x) and inverse it

    constexpr F a = 0.147;
    constexpr F inv_a = 1. / a;
    constexpr F tmp = (2 / (M_PI * a));
    const double ln1mx2 = Kokkos::log((1. - x) * (1. + x));
    const F term1 = tmp + (0.5 * ln1mx2);
    const F term2 = inv_a * ln1mx2;
    return Kokkos::copysign(Kokkos::sqrt(Kokkos::sqrt(term1 * term1 - term2) - term1), x);
  }

  // Inverse error function approximation (Using the rational approximation)
  //  template <FloatingPointType F> KOKKOS_INLINE_FUNCTION constexpr F erf_inv(F x)
  // {
  //     constexpr F a[4] = {0.147, 0.147, 0.147, 0.147};
  //     constexpr F b[4] = {-1.0, 0.5, -0.5, 1.0};
  //     KOKKOS_ASSERT(x <= -1.0 || x >= 1.0);

  //     F z = (x < 0.0) ? -x : x;

  //     F t = 2.0 / (Kokkos::numbers::pi * 0.147) + 0.5 * Kokkos::log(1.0 - z);
  //     F result = a[0] * Kokkos::pow(t, b[0]) +
  //                a[1] * Kokkos::pow(t, b[1]) +
  //                a[2] * Kokkos::pow(t, b[2]) +
  //                a[3] * Kokkos::pow(t, b[3]);
  //     KOKKOS_ASSERT(Kokkos::isfinite(result));
  //     return result;
  // }

  /**
  @brief Computes the inverse CDF (probit function) of a normal distribution.

  This function returns the quantile function of a normal distribution
  with mean `mean` and standard deviation `stddev`, given a probability `p`.

  @tparam F Floating-point type (must satisfy `FloatingPointType`).
  @param p Probability value in the range (0,1).
  @param mean Mean of the normal distribution.
  @param stddev Standard deviation of the normal distribution.
  @return The corresponding value `x` such that `P(X ≤ x) = p` for X ~ N(mean, stddev).

  @note This implementation uses the inverse error function `erfinv(x)`
        for accuracy and efficiency. Extreme values may be clamped for stability.

  @warning `p` must be strictly in (0,1), otherwise the result is undefined.
*/
  template <FloatingPointType F> KOKKOS_INLINE_FUNCTION F norminv(F p, F mean, F stddev)
  {
    constexpr F erfinv_lb = -5;
    constexpr F erfinv_up = 5;
    auto clamped_inverse = Kokkos::clamp(erfinv(2 * p - 1), erfinv_lb, erfinv_up);
    KOKKOS_ASSERT(Kokkos::isfinite(stddev * Kokkos::numbers::sqrt2 * clamped_inverse));
    return mean + stddev * Kokkos::numbers::sqrt2 * clamped_inverse;
  }

  /**
  @brief Computes the standard normal probability density function (PDF).

  The standard normal PDF is given by:
  \f[
      \phi(x) = \frac{1}{\sqrt{2\pi}} e^{-x^2 / 2}
  \f]

  @tparam F Floating-point type (must satisfy `FloatingPointType`).
  @param x Input value.
  @return Value of the standard normal PDF at `x`.

  @note This function is numerically stable.
*/
  template <FloatingPointType F> KOKKOS_INLINE_FUNCTION F std_normal_pdf(F x)
  {
    // NOLINTBEGIN(cppcoreguidelines-avoid-magic-numbers,readability-magic-numbers)
    constexpr double inv_sqrt_2_pi = 0.3989422804014327; // 1/sqrt(2pi)
    return inv_sqrt_2_pi * Kokkos::exp(-0.5 * x * x);
    // NOLINTEND(cppcoreguidelines-avoid-magic-numbers,readability-magic-numbers)
  }

  /**
  @brief Computes the standard normal cumulative distribution function (CDF).

  The standard normal CDF is given by:
  \f[
      \Phi(x) = \frac{1}{2} \left( 1 + \operatorname{erf} \left( \frac{x}{\sqrt{2}} \right) \right)
  \f]

  @tparam F Floating-point type (must satisfy `FloatingPointType`).
  @param x Input value.
  @return Value of the standard normal CDF at `x`, which is P(X ≤ x) for X ~ N(0,1).
*/
  template <FloatingPointType F> KOKKOS_INLINE_FUNCTION F std_normal_cdf(F x)
  {
    // NOLINTBEGIN(cppcoreguidelines-avoid-magic-numbers,readability-magic-numbers)
    return 0.5 * (1 + Kokkos::erf(x / Kokkos::numbers::sqrt2));
    // NOLINTEND(cppcoreguidelines-avoid-magic-numbers,readability-magic-numbers)
  }

  /*DISTRIBUTIONS*/

  /**
  @brief Represents a normal (Gaussian) probability distribution.

  The normal distribution is parameterized by a mean  and a standard deviation.
  It supports random sampling, computing statistical properties.

  @tparam F Floating-point type (must satisfy `FloatingPointType`).
*/
  template <FloatingPointType F> struct Normal
  {
    F mu;    ///< Mean
    F sigma; ///< Standard deviation

    /**
      @brief Draws a random sample from the distribution.
      @tparam DeviceType The Kokkos execution device.
      @param gen Kokkos random number generator.
      @return A normally distributed random sample.
    */
    template <class DeviceType>
    KOKKOS_INLINE_FUNCTION F draw(Kokkos::Random_XorShift1024<DeviceType>& gen) const
    {
      return draw_from(gen, mu, sigma);
    }

    /**
      @brief Static method to draw a sample from N(μ, σ).
      @tparam DeviceType The Kokkos execution device.
      @param gen Kokkos random number generator.
      @param mu Mean of the distribution.
      @param sigma Standard deviation of the distribution.
      @return A normally distributed random sample.
    */
    template <class DeviceType>
    static KOKKOS_INLINE_FUNCTION F draw_from(Kokkos::Random_XorShift1024<DeviceType>& gen,
                                              F mu,
                                              F sigma)
    {
      return gen.normal(mu, sigma);
    }

    /**
      @brief Returns the mean of the distribution.
      @return Mean (μ).
    */
    KOKKOS_INLINE_FUNCTION F mean() const
    {
      return mu;
    }

    /**
      @brief Returns the variance of the distribution.
      @return Variance (σ²).
    */
    KOKKOS_INLINE_FUNCTION F var() const
    {
      return sigma * sigma;
    }

    /**
      @brief Returns the standard deviation of the distribution.
      @return Standard deviation (σ).
    */
    KOKKOS_INLINE_FUNCTION F stddev() const
    {
      return sigma;
    }

    /**
      @brief Returns the skewness of the distribution.
      @return Skewness (always 0 for a normal distribution).
    */
    KOKKOS_INLINE_FUNCTION F skewness() const
    {
      return F(0);
    }
  };

  // To sample from [a, b] where |a - b| << 1, methods perform better
  // if we draw from xf = [a * factor, b * factor] and then scale back using x = xf / factor.
  template <FloatingPointType F> struct TruncatedNormal
  {

    F mu;    // Mean
    F sigma; // Standard deviation
    F lower; // Standard deviation
    F upper; // Standard deviation

    constexpr TruncatedNormal(F m, F s, F l, F u) : mu(m), sigma(s), lower(l), upper(u)
    {
      X_ASSERT(mu > lower);
      X_ASSERT(mu < upper);
      KOKKOS_ASSERT(mu > lower && mu < upper);
    }

    template <class DeviceType>
    KOKKOS_INLINE_FUNCTION F draw(Kokkos::Random_XorShift1024<DeviceType>& gen) const
    {
      return draw_from(gen, mu, sigma, lower, upper);
    }

    template <class DeviceType>
    static KOKKOS_INLINE_FUNCTION F
    draw_from(Kokkos::Random_XorShift1024<DeviceType>& gen, F mu, F sigma, F lower, F upper)
    {
      // NOLINTBEGIN(cppcoreguidelines-avoid-magic-numbers,readability-magic-numbers)
      const F rand = static_cast<F>(gen.drand());

      // Max bounded because if sigma <<1 z -> Inf not wanted because erf/erfc/erfinv are not stable
      // for extrem value Min bounded if |mu-bound| <<1 z -> 0 which is also not wanted for error
      // function

      F zl = Kokkos::clamp(lower - mu / sigma, F(1e-3), F(1e3));
      F zu = Kokkos::clamp(upper - mu / sigma, F(1e-3), F(1e3));

      F pl = 0.5 * Kokkos::erfc(-zl / Kokkos::numbers::sqrt2);
      KOKKOS_ASSERT(Kokkos::isfinite(pl)&&"Truncated normal draw leads is Nan of Inf with given parameters");
      F pu = 0.5 * Kokkos::erfc(-zu / Kokkos::numbers::sqrt2);
      KOKKOS_ASSERT(Kokkos::isfinite(pu)&&"Truncated normal draw leads is Nan of Inf with given parameters");
      F p = rand * (pu - pl) + pl;
      F x = norminv(p, mu, sigma);
      KOKKOS_ASSERT(Kokkos::isfinite(x)&&"Truncated normal draw leads is Nan of Inf with given parameters");
      return x;

      // NOLINTEND(cppcoreguidelines-avoid-magic-numbers,readability-magic-numbers)
    }

    KOKKOS_INLINE_FUNCTION F mean() const
    {
      const F alpha = (lower - mu) / sigma;
      const F beta = (upper - mu) / sigma;
      F Z = std_normal_cdf(beta) - std_normal_cdf(alpha);
      return mu + sigma * (std_normal_pdf(alpha) - std_normal_pdf(beta)) / Z;
    }

    KOKKOS_INLINE_FUNCTION F var() const
    {
      const F alpha = (lower - mu) / sigma;
      const F beta = (upper - mu) / sigma;
      F Z = std_normal_cdf(beta) - std_normal_cdf(alpha);
      const auto phi_b = std_normal_pdf(beta);
      const auto phi_a = std_normal_pdf(alpha);
      const auto tmp = (phi_a - phi_b) / Z;
      const auto tmp2 = (alpha * phi_a - beta * phi_b) / Z;
      return sigma * sigma * (1 - tmp2 - tmp * tmp);
    }
    KOKKOS_INLINE_FUNCTION F skewness() const
    {
      return 0.;
    }
  };

  template <FloatingPointType F> struct ScaledTruncatedNormal
  {

    F scale_factor;
    F inverse_factor;
    TruncatedNormal<F> dist;

    constexpr ScaledTruncatedNormal(F factor, F m, F s, F l, F u)
        : scale_factor(factor), inverse_factor(1. / scale_factor),
          dist(scale_factor * m, s * scale_factor, scale_factor * l, scale_factor * u)
    {
    }

    template <class DeviceType>
    static KOKKOS_INLINE_FUNCTION F draw_from(
        Kokkos::Random_XorShift1024<DeviceType>& gen, F factor, F mu, F sigma, F lower, F upper)
    {
      return 1. / factor *
             TruncatedNormal<F>::draw_from(
                 gen, factor * mu, factor * sigma, factor * lower, factor * upper);
    }

    template <class DeviceType>
    KOKKOS_INLINE_FUNCTION F draw(Kokkos::Random_XorShift1024<DeviceType>& gen) const
    {
      return inverse_factor * dist.draw(gen);
    }

    KOKKOS_INLINE_FUNCTION F mean() const
    {

      return inverse_factor * dist.mean();
    }

    KOKKOS_INLINE_FUNCTION F var() const
    {
      return (inverse_factor * inverse_factor) * dist.var();
    }
    KOKKOS_INLINE_FUNCTION F skewness() const
    {
      return 0.;
    }
  };

  template <FloatingPointType F> struct LogNormal
  {
    F mu;    // Mean
    F sigma; // Standard deviation

    template <class DeviceType>
    KOKKOS_INLINE_FUNCTION F draw(Kokkos::Random_XorShift1024<DeviceType>& gen) const
    {
      return Kokkos::exp(gen.normal(mu, sigma));
    }

    KOKKOS_INLINE_FUNCTION F mean() const
    {
      return Kokkos::exp(mu + sigma * sigma / 2);
    }
    KOKKOS_INLINE_FUNCTION F var() const
    {
      const auto sigma2 = sigma * sigma;
      return (Kokkos::exp(sigma2) - 1) * Kokkos::exp(2 * mu + sigma2);
    }
    KOKKOS_INLINE_FUNCTION F skewness() const
    {
      const auto sigma2 = sigma * sigma;
      return (Kokkos::exp(sigma2) + 2) * Kokkos::sqrt(Kokkos::exp(sigma2) - 1);
    }
  };

  template <FloatingPointType F> struct SkewNormal
  {
    F xi;    // Mean
    F omega; // Standard deviation
    F alpha;

    template <class DeviceType>
    KOKKOS_INLINE_FUNCTION F draw(Kokkos::Random_XorShift1024<DeviceType>& gen) const
    {

      const double Z0 = gen.normal();
      const double Z1 = gen.normal();
      const double delta = alpha / Kokkos::sqrt(1. + alpha * alpha);
      const double scale_factor = Kokkos::sqrt(1. + delta * delta);
      const double X = (Z0 + delta * Kokkos::abs(Z1)) / scale_factor;
      return xi + omega * X;
    }
    KOKKOS_INLINE_FUNCTION F mean() const
    {
      return xi + omega * (alpha / (Kokkos::sqrt(1 + alpha * alpha))) * Kokkos::sqrt(2 / M_PI);
    }
    KOKKOS_INLINE_FUNCTION F var() const
    {
      const auto delta = alpha / (Kokkos::sqrt(1 + alpha * alpha));
      return omega * omega * (1 - 2 * delta * delta / M_PI);
    }
    KOKKOS_INLINE_FUNCTION F skewness() const
    {
      const auto delta = alpha / (Kokkos::sqrt(1 + alpha * alpha));
      return ((4 - M_PI) * Kokkos::pow(delta * std::sqrt(2 / M_PI), 3)) /
             Kokkos::pow(1 - 2 * delta * delta / M_PI, 1.5);
    }
  };

} // namespace MC::Distributions

#endif