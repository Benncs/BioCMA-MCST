#include <cassert>
#include <common/traits.hpp>
#include <cstdint>

static_assert(FloatingPointType<double>, "double is float");
static_assert(FloatingPointType<const float&>, "cvref ");
static_assert(!FloatingPointType<int>, "int is not float");

static_assert(IntegerType<int>, "int is integer");
static_assert(IntegerType<const uint64_t&>, "cvref int");
static_assert(!IntegerType<bool>, "bool int");
static_assert(!IntegerType<double>, "double is not integer");
static_assert(!IntegerType<int*>, "pointer is not integer");

static_assert(NumberType<int>, "int is number");
static_assert(NumberType<float>, "float is number");
static_assert(!NumberType<bool>, "bool is not  number");

void
test_almost_equal()
{

  assert(almost_equal(1.0, 1.0));
  assert(almost_equal(1.0, 1.0 + 1e-20));
  assert(!almost_equal(1.0, 1.1));

  // Tolerance is inclusive
  assert(almost_equal(1.0, 1.05, 0.1));
  assert(!almost_equal(1.0, 1.05, 0.01));

  assert(almost_equal(-1.0, -1.0));
  assert(!almost_equal(-1.0, 1.0));

  // float
  assert(almost_equal(1.0F, 1.0F, 1e-6F));
  assert(!almost_equal(1.0F, 2.0F, 1e-6F));

  assert(!almost_equal(1., 2.0F, 1e-6F));

  // Mixed types compare in their common type, in both orders
  assert(almost_equal(1.0F, 1.0));
  assert(almost_equal(1.0, 1.0F));
  assert(!almost_equal(1.0F, 2.0));
  assert(!almost_equal(1.0, 2.0F));

  // Integer and floating point
  assert(almost_equal(1, 1.0, 1e-6));
  assert(!almost_equal(1, 2.0, 1e-6));

  // The tolerance is third type
  assert(almost_equal(1.0, 1.05F, 0.1));
  assert(!almost_equal(1.0, 1.05F, 0.01F));
}

void
test_almost_equal_ptr()
{
  const double a = 1.0;
  const double b = 1.0;
  const double c = 2.0;

  // Pointers to mixed types too
  const float f = 1.0F;

  assert(almost_equal(&a, &b));
  assert(!almost_equal(&a, &c));

  // Null on either side is false, never a dereference
  assert(!almost_equal(&a, static_cast<const double*>(nullptr)));
  assert(!almost_equal(static_cast<const double*>(nullptr), &a));
  assert(!almost_equal(static_cast<const double*>(nullptr),
                       static_cast<const double*>(nullptr)));

  assert(almost_equal(&a, &f));
  assert(!almost_equal(&c, &f));
}

int
main()
{
  test_almost_equal();
  test_almost_equal_ptr();
}
