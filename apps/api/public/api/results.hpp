#ifndef __API_RESULTS_HPP__
#define __API_RESULTS_HPP__

#include <string>
#include <variant>


struct Success
{
};
template <typename T> struct Result : protected std::variant<Success, T>
{
  explicit constexpr Result() noexcept : std::variant<Success, T>{Success{}} {};
  constexpr explicit Result(T const&& t) noexcept : std::variant<Success, T>{t}
  {
  }
  constexpr explicit operator bool() const noexcept
  {
    return valid();
  }

  [[nodiscard]] constexpr bool valid() const noexcept
  {
    return std::holds_alternative<Success>(*this);
  }
  [[nodiscard]] constexpr bool invalid() const noexcept
  {
    return !valid();
  }
  [[nodiscard]] constexpr auto get() const noexcept -> T
  {
    return (invalid() ? std::get<T>(*this) : T());
  }
};

struct ApiResult : Result<std::string>
{
  explicit ApiResult(std::string_view t) noexcept : Result<std::string>(std::string(t))
  {
  }

  explicit constexpr ApiResult() noexcept = default;
  constexpr int to_c_ret_code()
  {
    return (valid()) ? 0 : -1;
  }
};

#endif 