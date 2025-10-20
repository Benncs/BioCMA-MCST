#ifndef __MODEL_CONFIG_LOADER__
#define __MODEL_CONFIG_LOADER__

#include "Kokkos_Core.hpp"
#include "common/common.hpp"
#include "mc/prng/prng_extension.hpp"
#include <mc/traits.hpp>
#include <models/fixed_length.hpp>

namespace
{

}

namespace Models
{

  template <ConfigurableModel Model> Model::Config get_model_configuration()
  {
    return {}; // TODO
  }

  template <> inline FixedLength::Config get_model_configuration<FixedLength>()
  {
    using param_type = FixedLength::Params;

    float l_max = 2e-6;
    float l_men = 1e-6;
    float lvar = 2e-6 / 10.;
    char* ptr = std::getenv("LMAX");

    if (ptr != nullptr)
    {
      l_men = std::stof(ptr);
      std::cout << "Use env " << l_max << std::endl;
    }

    ptr = std::getenv("LVAR");
    if (ptr != nullptr)
    {
      lvar = std::stof(ptr);
      std::cout << "Use env " << lvar << std::endl;
    }

    auto l_initial_dist = MC::Distributions::Exponential<float>(l_men);

    auto host = Kokkos::View<param_type, HostSpace>("params");
    const auto params = param_type{l_initial_dist};
    host() = params;
    auto device = Kokkos::create_mirror_view_and_copy(
        ComputeSpace(), host, "params_device");
    return device;
  }

}; // namespace Models

#endif
