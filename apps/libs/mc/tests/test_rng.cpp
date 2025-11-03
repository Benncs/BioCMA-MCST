#include "common/traits.hpp"
#include <cassert>
#include <iostream>
#include <mc/prng/prng.hpp>

template <FloatingPointType T> void test_default_range(MC::KPRNG& rng)
{
  // Test 1: double_uniform() with default range
  auto result1 = rng.uniform<T>();
  assert(result1 > T(0.) && result1 < T(1.));
}

template <FloatingPointType T> void test_custom_range(MC::KPRNG& rng)
{
  // Test 2: double_uniform(a, b) with custom range
  T result = rng.uniform<T>(1.0, 3.0);
  assert(result > T(1.0) && result < T(3.0));
}

void test()
{
  MC::KPRNG rng;

  // Test 1: double_uniform() with default range
  test_default_range<float>(rng);
  test_default_range<double>(rng);

  // Test 2: double_uniform(a, b) with custom range
  test_custom_range<float>(rng);
  test_custom_range<double>(rng);

  // Test 3: dynamic size double_uniform(size_t n_sample, double a, double b)
  const size_t n_sample = 500;
  auto result3 = rng.double_uniform(n_sample, -1.0, 0.);

  for (size_t i = 0; i < n_sample; ++i)
  {
    assert(result3[i] > -1.0 && result3[i] < 0.);
  }

  // Test 4: double_uniform() with array size of n_r
  constexpr size_t n_sample_static = 25;
  const auto result4 = rng.double_uniform<n_sample_static>();
  for (auto&& i : result4)
  {
    assert(i > 0. && i < 1.);
  }

  const auto result5 = rng.random_view<double, n_sample_static>();
  for (size_t i = 0; i < n_sample_static; ++i)
  {
    std::cout << result5[i] << std::endl;
    assert(result5[i] > 0. && result5[i] < 1.);
  }

  // Test 5: uniform_u(uint64_t a, uint64_t b)
  uint64_t result6 = rng.uniform_u(10, 20);
  assert(result6 > 10 && result6 < 20);

  std::cout << "All tests passed!\n";
}

int main()
{

  Kokkos::initialize();
  test();
  Kokkos::finalize();
  return 0;
}
