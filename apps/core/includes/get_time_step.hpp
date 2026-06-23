#ifndef __CORE_GET_TIME_STEP_HPP__
#define __CORE_GET_TIME_STEP_HPP__

#include <cma_utils/alias.hpp>

double get_time_step(double user_deta_time,
                     const CmaUtils::TransitionnerPtrType& iterator);

#endif
