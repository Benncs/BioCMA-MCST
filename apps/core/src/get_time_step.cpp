#include <cma_utils/alias.hpp>
#include <get_time_step.hpp>
double
get_time_step(double user_deta_time,
              const CmaUtils::TransitionnerPtrType& iterator)
{

  // internal hydrodynamic time scales. To account for this, the simulation's
  // explicit time step is calculated to approximate a CFL condition, with the
  // formula: time_step = min(residence_time) / 100. This approach ensures
  // that the fluid movement between two steps is accurately represented
  // without losing flow information.
  double delta_time = user_deta_time;
  if (delta_time <= 0)
  {

    const auto min_residence_time = CmaUtils::get_min_residence_time(iterator);

    if (min_residence_time != std::numeric_limits<double>::max()
        && min_residence_time != 0.)
    {

      delta_time = min_residence_time / 100.;
    }
    else
    {
      // should throw  cause if delta <=0 we have incorrect timstep if
      // first branch fails delta_time is unchanged.
      // + If min_element fails flowmap might be invalid then default value
      // is not needed
      throw std::invalid_argument("No time step given and impossibe to "
                                  "estimate it with given flowmap");
    }
  }

  return delta_time;
}
