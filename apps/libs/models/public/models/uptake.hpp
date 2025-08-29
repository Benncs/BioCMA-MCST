#ifndef __MODELS_UPTAKES_HPP__
#define __MODELS_UPTAKES_HPP__

#include "common/traits.hpp"
#include <cstddef>
#include <mc/traits.hpp>

#define MONOD_RATIO(__c1__, __x1__, __k1__) ((__c1__) * (__x1__) / ((__x1__) + (__k1__)))

namespace Models
{


  template<FloatingPointType  F>
  struct UptakeDefault
  {
    using FloatType = F;
    // MODEL_CONSTANT FloatType NPermease_init = 1;
    // MODEL_CONSTANT FloatType NPermease_max = 200;
    MODEL_CONSTANT FloatType k = 1e-3;
    MODEL_CONSTANT FloatType beta = 5;
    MODEL_CONSTANT FloatType delta = 10;
    MODEL_CONSTANT FloatType tau_new_permease = 40.;
    MODEL_CONSTANT FloatType tau_rm_perm = 200.;
    MODEL_CONSTANT FloatType tau_pts = 20.;
    MODEL_CONSTANT FloatType tau_Au = 40.;
    MODEL_CONSTANT FloatType tau_Ad = 5.;
  };

  // CHECK_MODEL(Uptake<UptakeDefault<float>>)
} // namespace Models

#endif
