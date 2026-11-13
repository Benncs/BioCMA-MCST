#include <cassert>
#include <cmath>
#include <common/maths.hpp>
#include <common/traits.hpp>
#include <iostream>

// The fast path is a degree 3 fit on the mantissa, it is not accurate to
// float precision. These bounds are what it actually achieves, not a wish.
constexpr float abs_epsilon = 8e-3F;
constexpr float rel_epsilon = 1e-2F;

void
check(float x)
{
  const float fast = CommonMaths::_ln<false>(x);
  const float exact = CommonMaths::_ln<true>(x);
  const float absolute_error = std::fabs(fast - exact);
  const float relative_error = absolute_error / std::fabs(exact);

  if (absolute_error > abs_epsilon && relative_error > rel_epsilon)
  {
    std::cerr << "ln(" << x << ") fast=" << fast << " exact=" << exact
              << " abs=" << absolute_error << " rel=" << relative_error << '\n';
    assert(false);
  }
}

void
test_kokkos_path_is_std_log()
{
  // _ln<true> must be plain log, it is the reference the fast path is
  // measured against
  assert(almost_equal(CommonMaths::_ln<true>(2.0F), std::log(2.0F), 1e-7F));
  assert(almost_equal(CommonMaths::_ln<true>(0.5F), std::log(0.5F), 1e-7F));
}

void
test_agreement()
{
  // over the range the exponential draw uses, (0,1]
  for (int i = 1; i <= 10000; ++i)
  {
    check(static_cast<float>(i) / 10000.F);
  }

  check(1.F);
  check(2.F);
  check(10.F);
  check(1000.F);
  check(1e6F);
}

void
test_exact_at_one()
{
  assert(std::fabs(CommonMaths::_ln<false>(1.F)) < abs_epsilon);
  assert(CommonMaths::_ln<true>(1.F) == 0.F);
}

int
main()
{
  test_kokkos_path_is_std_log();
  test_agreement();
  test_exact_at_one();
}
