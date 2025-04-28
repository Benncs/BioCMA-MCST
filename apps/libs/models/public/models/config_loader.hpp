#ifndef __MODEL_CONFIG_LOADER__
#define __MODEL_CONFIG_LOADER__

#include <mc/traits.hpp>


namespace Models
{

  template <ConfigurableModel Model> Model::Config get_model_configuration()
  {
    return {}; // TODO
  }
}; // namespace Models

#endif