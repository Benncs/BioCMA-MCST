#ifndef __API_RESULTS_HPP__
#define __API_RESULTS_HPP__

#include <common/results.hpp>
#include <string>

struct ApiResult : Result<Success, std::string>
{
  explicit ApiResult(std::string_view t) noexcept
      : Result<Success, std::string>(std::string(t))
  {
  }

  explicit constexpr ApiResult() noexcept = default;
  constexpr int to_c_ret_code()
  {
    return (valid()) ? 0 : -1;
  }
};

#endif