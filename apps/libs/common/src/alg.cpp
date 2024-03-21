#include <cmath>
#include <common/alg.hpp>

double naive_newton(const std::function<double(double)> &f,
                    double x_init,
                    bool *success,
                    double tolerance)
{
  double x = x_init;
  constexpr double newton_dx = 1e-8;
  auto d = [&](double eps_guess) -> double
  {
    return (f(eps_guess + newton_dx) - f(eps_guess - newton_dx)) /
           (2 * newton_dx);
  };

  for (auto iter = 0; iter < 50; ++iter)
  {
    double val = f(x);
    if (std::isnan(val))
    {
      if (success)
      {
        *success = false;
        return x;
      }
    }
    
    double derivative_val = d(x);

    double y = val / derivative_val;
    x -= y;

    if (std::abs(val) < tolerance)
    {

      if (success)
        *success = true;
      return x;
    }
  }

  if (success)
    *success = false;
  return x;
}