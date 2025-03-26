#ifndef __COMMON_RESULTS_HPP__
#define __COMMON_RESULTS_HPP__
#include <stdexcept>
#include <utility>
#include <variant>

struct Success
{
};
// TODO WIP
template <typename S, typename T> struct Result : protected std::variant<S, T>
{
  explicit constexpr Result() noexcept : std::variant<S, T>{S{}} {};
  constexpr explicit Result(T const&& t) noexcept : std::variant<S, T>{t}
  {
  }

  constexpr explicit Result(S&& s) noexcept : std::variant<S, T>{std::move(s)}
  {
  }

  constexpr explicit operator bool() const noexcept
  {
    return valid();
  }

  [[nodiscard]] constexpr bool valid() const noexcept
  {
    return std::holds_alternative<S>(*this);
  }
  [[nodiscard]] constexpr bool invalid() const noexcept
  {
    return !valid();
  }

  [[nodiscard]] constexpr auto get() const noexcept -> T
  {
    return (invalid() ? std::get<T>(*this) : T());
  }

  [[nodiscard]] auto gets() const -> S
  {
    if (!invalid())
    {
      return std::move(std::get<S>(*this));
    }
    throw std::runtime_error("Deref None");
  }

  // template <typename Func, typename Err> auto match_const(Func&& f, Err&& r) const noexcept
  // {
  //   if (this->valid())
  //   {
  //     return f(std::get<S>(*this));
  //   }
  //   else
  //   {
  //     return r(std::get<T>(*this));
  //   }
  // }

  template <typename Func, typename Err> auto match(Func&& f, Err&& r) noexcept
  {
    if (this->valid())
    {
      return f(std::forward<S>(std::get<S>(*this)));
    }

    return r(std::forward<T>(std::get<T>(*this)));
  }
};

#endif