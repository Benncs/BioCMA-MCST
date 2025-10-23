
#include "Kokkos_Core.hpp"
#include <iostream>
#include <random>
#include <simulation/probability_leaving.hpp>
#ifdef NDEBUG
#  undef NDEBUG
#endif
#include <cassert>
float random_float(float min, float max)
{
  static std::random_device rd;
  static std::mt19937 gen(rd());
  std::uniform_real_distribution<float> dis(min, max);
  return dis(gen);
}

double log_fast_ankerl(double a)
{
  static_assert(__BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__,
                "Little endian is required!");
  union
  {
    double d;
    int x[2];
  } u = {a};
  return (u.x[1] - 1072632447) * 6.610368362777016e-7;
}

void test_ln()
{
  float abs_epsilon = 1e-3; // Allowable error
  float rel_epsilon = 0.7;  // Allowable error

  for (int i = 0; i < 1e6; ++i)
  {
    auto x = random_float(0, 1);
    float computed_no_kokkos = Simulation::KernelInline::_ln<false>(x);
    float computed_with_kokkos = Simulation::KernelInline::_ln<true>(x);
    std::cout << computed_no_kokkos << " " << computed_with_kokkos << std::endl;
    float absolute_error = std::fabs(computed_no_kokkos - computed_with_kokkos);
    float relative_error = (absolute_error / std::fabs(computed_with_kokkos));
    assert(absolute_error < abs_epsilon);
    assert(relative_error < rel_epsilon);
  }
}

int main()
{
  Kokkos::initialize();
  test_ln();
  Kokkos::finalize();
}