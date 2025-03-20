#ifndef __PMESSAGE_HPP__
#define __PMESSAGE_HPP__

#include <type_traits>
namespace WrapMPI
{
  // static bool is_mpi_init = false;
  int critical_error() noexcept;
  void barrier() noexcept;
  void finalize() noexcept;
  bool is_initialized() noexcept;
  enum class SIGNALS : char
  {
    STOP,
    RUN,
    NOP,
    DUMP
  };

    template <typename T>
  concept POD_t =
      std::is_standard_layout_v<T> && std::is_trivially_copyable_v<T> &&
      std::is_trivially_destructible_v<T> && std::is_trivially_default_constructible_v<T>;


} // namespace WrapMPI

#endif //__PMESSAGE_HPP__
