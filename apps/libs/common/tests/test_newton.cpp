#include <cassert>
#include <cmath>
#include <common/common.hpp>
#include <iostream>

void validation_test()
{
  auto sqrt_function = [](double x) { return x * x - 4; };
  bool success = false;
  double result = naive_newton(sqrt_function, 19.0, &success);
  std::cerr << success << " " << result << '\n';
  if (success)
  {
    assert(std::abs(result - 2.0) < 1e-6);
  }
  else
  {
    assert(false);
  }
}

void validation_test_2()
{
  auto invalid_function = [](double x) { return std::log(x); };
  bool success = false;
  double result = naive_newton(invalid_function, 0.0, &success);
  std::cerr << success << " " << result << '\n';
  assert(!success && std::isnan(result));
}

int main()
{
  validation_test();
  validation_test_2();
}