#ifndef __CMA_UTILS_HPP__
#define __CMA_UTILS_HPP__

#include <cma_utils/cache_hydro_state.hpp>

namespace CmaUtils
{
  struct IterationState
  {
    CmaUtils::PreCalculatedHydroState* liq;
    CmaUtils::PreCalculatedHydroState* gas;
  };
} // namespace CmaUtils

#endif