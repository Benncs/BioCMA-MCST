#include "Kokkos_Macros.hpp"
#include "Kokkos_MathematicalFunctions.hpp"
#include <models/config_loader.hpp>

#include <models/fixed_length.hpp>

#include <Kokkos_sampling/metropolis.hpp>

namespace Models
{
  FixedLength::Config
  FixedLength::get_config(const std::size_t n)
  {
    using float_type = Self::FloatType;
    Kokkos::View<float_type*, ComputeSpace> samples("samples", n);

    float_type lambda = Kokkos::log(2.F) / 1e-6F;

    lambda = Common::read_env_or("VLAMBDA", lambda);

    const std::size_t nbin = Common::read_env_or("VNBIN", 50);

    if (lambda != 0.F)
    {
      // Given that hist(samples)should be : lambda * Kokkos::exp(-lambda * x);
      // Use inverse sampling strategy within bin

      // Expected extra of the given distribution
      float_type e_min = Kokkos::exp(-lambda * Self::l_min_m);
      float_type e_max = Kokkos::exp(-lambda * Self::l_max_m);

      // Inverse exponential sapling between emin/emax
      // Same strategy as inverse sampling with u within [0,1]
      // but here u is not random,
      auto target = KOKKOS_LAMBDA(const FixedLength::FloatType y)
      {
        return -(1.0F / lambda) * Kokkos::log(e_min - y * (e_min - e_max));
      };

      // Prepare n_blocks so that n_blocks*nbin=n

      std::size_t n_blocks = (n + nbin - 1) / nbin;
      Kokkos::RangePolicy<ComputeSpace> pl0(0, n_blocks);
      // iterate of block and span all bins in each iteration
      // Then automatically place particle in the correct bin
      Kokkos::parallel_for(
          "init_fixed_lenght", pl0, KOKKOS_LAMBDA(const int i) {
            std::size_t base = i * nbin;
            for (std::size_t j = 0; j < nbin; ++j)
            {
              std::size_t index = base + j;
              if (index < n)
              {
                // Need to retrieve the bin position for the index
                // Inverse samplijng requires index in [0,1]
                // take density at the middle of the bin
                float_type q = (static_cast<float_type>(index) + 0.5F)
                               / static_cast<float_type>(n);

                samples(index) = target(q);
              }
            }
          });
    }
    else
    {
      throw std::runtime_error("Unimplemented yet");
    }

    return samples;
  }

  // FixedLength::Config
  // FixedLength::get_config(const std::size_t n)
  // {
  //   using float_type = Self::FloatType;
  //   float_type lambda = Kokkos::log(2.F) / 1e-6F;
  //   char* lambda_env = std::getenv("VLAMBDA");

  //   if (lambda_env != nullptr)
  //   {
  //     lambda = static_cast<float_type>(std::stod(lambda_env));
  //     Kokkos::printf("[Config] use env value %f\r\n", lambda);
  //   }

  //   int rc = 0;
  //   Kokkos::View<float_type*, ComputeSpace> samples("samples", n);
  //   if (lambda != 0.)
  //   {
  //     auto target = KOKKOS_LAMBDA(const FixedLength::FloatType x)
  //     {
  //       return lambda * Kokkos::exp(-lambda * x);
  //     };

  //     rc = Sampling::metropolis(target, samples, Self::l_min_m,
  //     Self::l_max_m);
  //   }
  //   else
  //   {
  //     constexpr FixedLength::FloatType mu
  //         = 1.3e-6; // maxlength is 2 and min i 1 so let be between
  //     constexpr FixedLength::FloatType sigma = 0.1;

  //     auto target_distribution = KOKKOS_LAMBDA(const FixedLength::FloatType
  //     x)
  //     {
  //       if (x <= 0.0)
  //       {
  //         return 0.0;
  //       }
  //       auto coefficient = 1.0 / (x * sigma * std::sqrt(2.0 * M_PI));
  //       auto exponent = -std::pow(std::log(x) - mu, 2) / (2.0 * sigma *
  //       sigma);

  //       return coefficient * std::exp(exponent);
  //     };

  //     rc = Sampling::metropolis(
  //         target_distribution, samples, Self::l_min_m, Self::l_max_m);
  //   }

  //   if (rc != 0)
  //   {
  //     throw std::runtime_error("FixedLength init: Error when sampling");
  //   }

  //   return samples;
  // }

} // namespace Models
