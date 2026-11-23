#include <Kokkos_Array.hpp>
#include <cassert>
#include <cmath>
#include <common/maths.hpp>
#include <common/static_for.hpp>
#include <common/traits.hpp>
#include <cstddef>

void
test_static_for_order()
{
  constexpr std::size_t n = 5;
  Kokkos::Array<std::size_t, n> seen{};
  std::size_t calls = 0;

  Common::static_for<n>(
      [&](auto i)
      {
        seen[calls] = decltype(i)::value;
        ++calls;
      });

  assert(calls == n);
  for (std::size_t k = 0; k < n; ++k)
  {
    assert(seen[k] == k);
  }
}

void
test_static_for_same_as_loop()
{
  constexpr std::size_t n = 7;
  Kokkos::Array<double, n> x{};
  for (std::size_t k = 0; k < n; ++k)
  {
    x[k] = 0.1 * static_cast<double>(k + 1);
  }

  double loop_sum = 0.;
  for (std::size_t k = 0; k < n; ++k)
  {
    loop_sum += x[k] * x[k];
  }

  double static_sum = 0.;
  Common::static_for<n>(
      [&](auto i)
      {
        constexpr std::size_t k = decltype(i)::value;
        static_sum += x[k] * x[k];
      });

  assert(static_sum == loop_sum);
}

void
test_static_for_index_is_constant()
{
  constexpr std::size_t n = 6;
  std::size_t n_even = 0;

  Common::static_for<n>(
      [&](auto i)
      {
        constexpr std::size_t k = decltype(i)::value;
        static_assert(k < n, "index out of the static_for range");
        if constexpr (k % 2 == 0)
        {
          ++n_even;
        }
      });

  assert(n_even == 3);
}

void
test_static_for_constexpr()
{
  constexpr std::size_t sum = []
  {
    std::size_t s = 0;
    Common::static_for<4>([&](auto i) { s += decltype(i)::value; });
    return s;
  }();

  static_assert(sum == 6, "static_for must run in a constant expression");
}

void
test_static_for_empty()
{
  std::size_t calls = 0;
  Common::static_for<0>([&](auto /*i*/) { ++calls; });
  assert(calls == 0);
}

int
main()
{
  test_static_for_order();
  test_static_for_same_as_loop();
  test_static_for_index_is_constant();
  test_static_for_constexpr();
  test_static_for_empty();
}
