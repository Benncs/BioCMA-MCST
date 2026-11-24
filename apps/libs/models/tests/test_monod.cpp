#include <cassert>
#include <common/traits.hpp>
#include <concepts>
#include <models/utils.hpp>

template <typename T>
  requires std::floating_point<T>
void
test_monod_batch_three()
{
  const T c_s = 1e-2F;
  const T k_s = 1e-3F;
  const T c_o2 = 5e-5F;
  const T k_o = 1e-5F;
  const T c_ac = 4e-4F;
  const T k_a = 1e-4F;

  const auto [m_s, m_o2, m_ac]
      = Models::monod_batch(c_s, k_s, c_o2, k_o, c_ac, k_a);
  T tol = T{ 1e-6 };
  assert(almost_equal(m_s, c_s / (k_s + c_s), tol));
  assert(almost_equal(m_o2, c_o2 / (k_o + c_o2), tol));
  assert(almost_equal(m_ac, c_ac / (k_a + c_ac), tol));
}

template <typename T>
  requires std::floating_point<T>
void
test_monod_batch_one()
{
  const T a = 3.;
  const T ka = 1.;

  const auto [m_a] = Models::monod_batch(a, ka);
  T tol = T{ 1e-15 };
  assert(almost_equal(m_a, a / (ka + a), tol));
}

int
main()
{
  test_monod_batch_one<float>();
  test_monod_batch_one<double>();
  test_monod_batch_three<float>();
  test_monod_batch_three<double>();
}
