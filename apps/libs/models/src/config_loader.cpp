#include <models/config_loader.hpp>

#include <models/fixed_length.hpp>

#include <Kokkos_sampling/metropolis.hpp>

namespace Models
{

  FixedLength::Config FixedLength::get_config(const std::size_t n)
  {

    using float_t = Self::FloatType;

    float_t lambda = 2e6;
    char* lambda_env = std::getenv("VLAMBDA");
    if (lambda_env != nullptr)
    {
      lambda = static_cast<float_t>(std::stod(lambda_env));
      Kokkos::printf("[Config] use env value %f\r\n", lambda);
    }

    auto target = KOKKOS_LAMBDA(const FixedLength::FloatType x)
    {
      return lambda * Kokkos::exp(-lambda * x);
    };

    Kokkos::View<float_t*, ComputeSpace> samples("samples", n);

    auto rc =
        Sampling::metropolis(target, samples, Self::l_min_m, Self::l_max_m);

    if (rc != 0)
    {
      throw std::runtime_error("FixedLength init: Error when sampling");
    }

    return samples;
  }

} // namespace Models
