#ifndef __COMMON_TRAITS_HPP__
#define __COMMON_TRAITS_HPP__
#include <cassert>
#include <cmath>
#include <type_traits>

#if defined __GNUC__
#  define LIKELY(EXPR) __builtin_expect(!!(EXPR), 1)
#else
#  define LIKELY(EXPR) (!!(EXPR))
#endif

#if defined NDEBUG
#  define X_ASSERT(CHECK) void(0)
#else
#  define X_ASSERT(CHECK)                                                      \
    (LIKELY(CHECK) ? void(0) : [] { assert(!(#CHECK)); }())
#endif

template <typename T>
concept FloatingPointType = std::is_floating_point_v<std::remove_cvref_t<T>>;

static_assert(FloatingPointType<double>, "double ok");
static_assert(FloatingPointType<float>, "float ok");
static_assert(!FloatingPointType<int>, "int ok");

constexpr double tolerance_equality_float = 1e-15;

template <typename T>
concept IntegerType = requires(T n) {
  requires std::is_integral_v<std::remove_cvref_t<T>>;
  requires !std::is_same_v<std::remove_cvref_t<T>, bool>;
  requires std::is_arithmetic_v<decltype(n + 1)>;
  requires !std::is_pointer_v<std::remove_cvref_t<T>>;
};

template <typename T>
concept NumberType = IntegerType<T> || FloatingPointType<T>;

// One parameter per argument, so mixed types compare in their common type.
// Numbers are cheap to copy, by value covers lvalue and rvalue alike
template <NumberType T, NumberType U, NumberType Tol = double>
inline bool
almost_equal(T val, U val2, Tol tolerance = tolerance_equality_float)
{
  using CommonT = std::common_type_t<T, U>;
  return std::abs(static_cast<CommonT>(val) - static_cast<CommonT>(val2))
         < static_cast<CommonT>(tolerance);
}

// Overload for pointers
template <NumberType T, NumberType U, NumberType Tol = double>
inline bool
almost_equal(const T* val,
             const U* val2,
             Tol tolerance = tolerance_equality_float)
{
  if (val == nullptr || val2 == nullptr)
  {
    return false; // Null pointer check
  }
  return almost_equal(*val, *val2, tolerance);
}

#endif
