#ifndef __COMMON_EXPORT_HPP__
#define __COMMON_EXPORT_HPP__

#include <common/alg.hpp>
#include <common/common_types.hpp>
#include <common/execinfo.hpp>
#include <common/simulation_parameters.hpp>
#include <biocma_cst_config.hpp>

#ifdef ENABLE_KOKKOS_PROFILING
#include <Kokkos_Profiling_ScopedRegion.hpp>
#define PROFILE_SECTION(__label_section__) Kokkos::Profiling::ScopedRegion region(__label_section__);
#else 
#define PROFILE_SECTION(__label_section__)  ;
#endif 

#ifndef NDEBUG
#define WARN_EXPERIMENTAL _Pragma("message( __FILE__ \" contains experimental Kokkos feature\")") // NOLINT
#else
#define WARN_EXPERIMENTAL

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
concept NumberType =
    requires(T n) { requires IntegerType<T> || std::is_floating_point_v<T>; };


template<NumberType T>
inline bool almost_equal(T val,T val2,T tolerance = tolerance_equality_float)
{
  return std::abs(val-val2)<tolerance;
}

#define EIGEN_INDEX(__VALUE__) static_cast<int>(__VALUE__)


#endif //__COMMON_EXPORT_HPP__