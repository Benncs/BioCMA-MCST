#ifndef __MODEL_CONFIG_LOADER__
#define __MODEL_CONFIG_LOADER__

#include "Kokkos_Core.hpp"
#include "common/common.hpp"
#include "mc/prng/prng_extension.hpp"
#include <kokkos-sampling/metropolis.hpp>
#include <mc/traits.hpp>
#include <models/fixed_length.hpp>
#include <stdexcept>

namespace Models
{

  template <ConfigurableModel Model>
  Model::Config get_model_configuration(std::size_t n)
  {
    return {}; // TODO
  }

  template <>
  inline FixedLength::Config get_model_configuration<FixedLength>(std::size_t n)
  {
    using float_t = FixedLength::FloatType;

    float_t lambda = 0;
    char* lambda_env = std::getenv("VLAMBDA");
    if (lambda_env != nullptr)
    {
      lambda = static_cast<float_t>(std::stod(lambda_env));
    }

    auto target = KOKKOS_LAMBDA(const FixedLength::FloatType x)
    {
      return lambda * Kokkos::exp(-lambda * x);
    };
    Kokkos::View<float_t*, ComputeSpace> samples("samples", n);

    auto rc = metropolis(target, samples, float_t(1e-6), float_t(2e-6));

    if (rc != 0)
    {
      throw std::runtime_error("FixedLength init: Error when sampling");
    }

    return samples;
  }

}; // namespace Models

#endif
