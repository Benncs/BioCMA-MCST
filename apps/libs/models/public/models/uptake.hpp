#ifndef __MODELS_UPTAKES_HPP__
#define __MODELS_UPTAKES_HPP__

#include "common/traits.hpp"
#include <cstddef>
#include <mc/traits.hpp>

#define MONOD_RATIO(__c1__, __x1__, __k1__)                                    \
  ((__c1__) * (__x1__) / ((__x1__) + (__k1__)))

namespace Models
{

  template <FloatingPointType F> struct UptakeDefault
  {
    using FloatType = F;
    // MODEL_CONSTANT FloatType NPermease_init = 1;
    // MODEL_CONSTANT FloatType NPermease_max = 200;
    MODEL_CONSTANT FloatType tau_ap_1 = 500;
    MODEL_CONSTANT FloatType tau_ap_2 = 200;
    MODEL_CONSTANT FloatType tau_ap_3 = 1000;

    MODEL_CONSTANT FloatType k = 1e-2;
    MODEL_CONSTANT FloatType beta = 5;
    MODEL_CONSTANT FloatType k_perm = 1e-3; // TODO
  };

  // CHECK_MODEL(Uptake<UptakeDefault<float>>)
} // namespace Models

#endif
