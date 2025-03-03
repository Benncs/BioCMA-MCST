#ifndef __PRNG_EXTENSION_HPP__
#define __PRNG_EXTENSION_HPP__

#include "Kokkos_Assert.hpp"
#include "Kokkos_MathematicalConstants.hpp"
#include <Kokkos_Random.hpp>
#include <cmath>
#include <common/traits.hpp>

namespace MC::Distributions
{

  template <typename T, typename F, class DeviceType>
  concept ProbabilityLaw =
      FloatingPointType<F> && requires(const T& obj, Kokkos::Random_XorShift64<DeviceType>& gen) {
        { obj.draw(gen) } -> std::same_as<F>;
        { obj.mean() } -> std::same_as<F>;
        { obj.var() } -> std::same_as<F>;
        { obj.skewness() } -> std::same_as<F>;
      };

  // Approximate inverse error function (needed for inverse CDF sampling)
  template <FloatingPointType F> KOKKOS_INLINE_FUNCTION F erfinv(F x)
  {

    // Approximation by Abramowitz and Stegun handbook, not really optimised
    // Check this for GPU impl https://people.maths.ox.ac.uk/gilesm/codes/erfinv/gems.pdf

    constexpr F a = 0.147;
    constexpr F inv_a = 1. / a;
    constexpr F tmp = (2 / (M_PI * a));
    // const double ln1mx2 = Kokkos::log((x-1.)*(1.+x));
    const F ln1mx2 = Kokkos::log(1. - x * x);
    const F term1 = tmp + (0.5 * ln1mx2);
    const F term2 = inv_a * ln1mx2;
    return Kokkos::copysign(Kokkos::sqrt(Kokkos::sqrt(term1 * term1 - term2) - term1), x);
  }

  // Inverse CDF (Probit function)
  template <FloatingPointType F> KOKKOS_INLINE_FUNCTION F norminv(F p, F mean, F stddev)
  {
    return mean + stddev * Kokkos::numbers::sqrt2 * erfinv(2 * p - 1);
  }

  template <FloatingPointType F> KOKKOS_INLINE_FUNCTION F std_normal_pdf(F x)
  {
    constexpr double inv_sqrt_2_pi = 0.3989422804014327; // 1/sqrt(2pi)

    return inv_sqrt_2_pi * Kokkos::exp(-0.5 * x * x);
  }

  template <FloatingPointType F> KOKKOS_INLINE_FUNCTION F std_normal_cdf(F x)
  {
    return 0.5 * (1 + Kokkos::erf(x / Kokkos::numbers::sqrt2));
  }

  template <FloatingPointType F> struct Normal
  {
    F mu;    // Mean
    F sigma; // Standard deviation

    template <class DeviceType>
    KOKKOS_INLINE_FUNCTION F draw(Kokkos::Random_XorShift64<DeviceType>& gen) const
    {
      return draw_from(gen, mu, sigma);
    }

    template <class DeviceType>
    static KOKKOS_INLINE_FUNCTION F draw_from(Kokkos::Random_XorShift64<DeviceType>& gen,
                                              F mu,
                                              F sigma)
    {
      return gen.normal(mu, sigma);
    }

    KOKKOS_INLINE_FUNCTION F mean() const
    {
      return mu;
    }
    KOKKOS_INLINE_FUNCTION F var() const
    {
      return sigma * sigma;
    }
    KOKKOS_INLINE_FUNCTION F skewness() const
    {
      return 0.;
    }
  };

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
    KOKKOS_INLINE_FUNCTION F draw(Kokkos::Random_XorShift64<DeviceType>& gen) const
    {
      return draw_from(gen, mu, sigma, lower, upper);
    }

    template <class DeviceType>
    static KOKKOS_INLINE_FUNCTION F
    draw_from(Kokkos::Random_XorShift64<DeviceType>& gen, F mu, F sigma, F lower, F upper)
    {
      const F rand = static_cast<F>(gen.frand());
      F pl = F(0.5) * (1 + Kokkos::erf((lower - mu) / (sigma * Kokkos::numbers::sqrt2)));
      F pu = F(0.5) * (1 + Kokkos::erf((upper - mu) / (sigma * Kokkos::numbers::sqrt2)));
      F p = rand * (pu - pl) + pl;
      return norminv(p, mu, sigma);
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

  template <FloatingPointType F> struct LogNormal
  {
    F mu;    // Mean
    F sigma; // Standard deviation

    template <class DeviceType>
    KOKKOS_INLINE_FUNCTION F draw(Kokkos::Random_XorShift64<DeviceType>& gen) const
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
    KOKKOS_INLINE_FUNCTION F draw(Kokkos::Random_XorShift64<DeviceType>& gen) const
    {
      //   double alpha2 =alpha * alpha;
      //   double omega2 = omega*omega;
      //   double m = Kokkos::sqrt((1+alpha2)/2);
      //   double U1 = gen.normal(m,omega2);
      //   double U2 = gen.normal(m,omega2);
      //   double U = Kokkos::max(U1,U2);
      //   double V = Kokkos::min(U1,U2);
      //   double lambda1 = (1+alpha) / Kokkos::sqrt(1 + alpha2);
      //   double lambda2 = (1-alpha) / Kokkos::sqrt(1  +alpha2);

      //   return U*lambda1+V*lambda2;

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