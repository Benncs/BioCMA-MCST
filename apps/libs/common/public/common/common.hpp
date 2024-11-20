#ifndef __COMMON_EXPORT_HPP__
#define __COMMON_EXPORT_HPP__

#include <biocma_cst_config.hpp>
#include <common/alg.hpp>
#include <common/common_types.hpp>
#include <common/execinfo.hpp>
#include <string_view>


#ifdef ENABLE_KOKKOS_PROFILING
#  include <Kokkos_Profiling_ScopedRegion.hpp>
#  define PROFILE_SECTION(__label_section__) Kokkos::Profiling::ScopedRegion region(__label_section__);
#else
#  define PROFILE_SECTION(__label_section__) ;
#endif

#ifndef NDEBUG
#  define WARN_EXPERIMENTAL _Pragma("message( __FILE__ \" contains experimental Kokkos feature\")") // NOLINT
#else
#  define WARN_EXPERIMENTAL

#endif
#define EIGEN_INDEX(__VALUE__) static_cast<int>(__VALUE__)

#ifndef NDEBUG
#  include <iostream>
#  include <source_location>
class Canary
{
public:
  Canary() = delete;
  Canary(const Canary &) = delete;
  Canary(Canary &&) = delete;
  Canary &operator=(const Canary &) = delete;
  Canary &operator=(Canary &&) = delete;

  explicit Canary(std::string_view lbl, std::source_location location) : _lbl(lbl)
  {
    // std::cout << "\033[1;31m" << location.function_name() << "\033[0m " << lbl << std::endl;
    std::cout << location.function_name()<<": " << lbl << std::endl;
  }

  ~Canary()
  {
    std::cout << "END " << _lbl << std::endl;
  }

private:
  std::string _lbl;
};
#  define MkCanary(x) Canary(x, std::source_location());
#else
#  define WARN_EXPERIMENTAL
#  define MkCanary(x)
#endif

void test_common();

double constexpr tolerance_equality_float = 1e-15;
template <typename T>
concept IntegerType = requires(T n) {
  requires std::is_integral_v<T>;
  requires !std::is_same_v<bool, T>;
  requires std::is_arithmetic_v<decltype(n + 1)>;
  requires !std::is_pointer_v<T>;
};

template <typename T>
concept NumberType = requires(T n) { requires IntegerType<T> || std::is_floating_point_v<T>; };

template <NumberType T> inline bool almost_equal(T val, T val2, T tolerance = tolerance_equality_float)
{
  return std::abs(val - val2) < tolerance;
}



#endif //__COMMON_EXPORT_HPP__