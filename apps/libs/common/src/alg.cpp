#include <cmath>
#include <common/alg.hpp>

#include <cmath>
#include <functional>

double naive_newton(const std::function<double(double)> &f,
                    double x_init,
                    bool *success,
                    double tolerance)
{
  constexpr int max_iterations = 50;
  constexpr double newton_dx = 1e-8;

  auto d = [&](double eps_guess) -> double
  {
    return (f(eps_guess + newton_dx) - f(eps_guess - newton_dx)) /
           (2 * newton_dx);
  };

  for (int iter = 0; iter < max_iterations; ++iter)
  {
    double val = f(x_init);
    if (std::isnan(val))
    {
      if (success)
        *success = false;
      return x_init;
    }

    double derivative_val = d(x_init);
    if (std::abs(derivative_val) < std::numeric_limits<double>::epsilon())
    {
      if (success)
        *success = false;
      return x_init;
    }

    double y = val / derivative_val;
    x_init -= y;

    if (std::abs(val) < tolerance)
    {
      if (success)
        *success = true;
      return x_init;
    }
  }

  if (success)
    *success = false;
  return x_init;
}
