#ifndef __COMMON_STATIC_FOR_HPP__
#define __COMMON_STATIC_FOR_HPP__

#include <Kokkos_Macros.hpp>
#include <cstddef>
#include <type_traits>
#include <utility>

namespace Common
{
  namespace impl
  {
    template <typename F, std::size_t... I>
    KOKKOS_FORCEINLINE_FUNCTION constexpr void
    static_for_impl(F& f, std::index_sequence<I...> /*unused*/)
    {
      (f(std::integral_constant<std::size_t, I>{}), ...);
    }
  } // namespace impl

  /**
   * @brief Compile time loop: calls f(integral_constant<0>) ...
   * f(integral_constant<N-1>), in order.
   *
   * Unrolled by construction (fold expression), so it does not rely on the
   * optimiser, on the host as on the device. The index reaches f as an
   * integral_constant: inside f, `decltype(i)::value` is a constant
   * expression (usable in `if constexpr` or as a template argument).
   */
  template <std::size_t N, typename F>
  KOKKOS_FORCEINLINE_FUNCTION constexpr void
  static_for(F&& f)
  {
    impl::static_for_impl(f, std::make_index_sequence<N>{});
  }
} // namespace Common

#endif
