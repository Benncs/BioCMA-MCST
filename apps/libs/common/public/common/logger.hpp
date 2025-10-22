#ifndef __COMMON_LOGGER_HPP__
#define __COMMON_LOGGER_HPP__

#include <source_location>
#include <string>
#include <string_view>

namespace IO
{

  namespace AnsiCode
  {
    constexpr auto blue = "\033[34m";
    constexpr auto red = "\033[31m";
    constexpr auto reset = "\033[0m";
    constexpr auto green = "\033[32m";
  } // namespace AnsiCode

  namespace Unicode
  {
    constexpr auto red_circle = "\U0001F534";
  } // namespace Unicode

  class Logger
  {
  public:
    Logger() = default;

    Logger(const Logger&) = default;
    Logger(Logger&&) = default;

    Logger& operator=(const Logger&) = default;
    Logger& operator=(Logger&&) = default;

    virtual ~Logger() = default;

    virtual void debug(std::string_view message) = 0;

    virtual void print(std::string_view prefix, std::string_view message) = 0;

    virtual void alert(std::string_view prefix, std::string_view message) = 0;

    virtual void
    error(std::string_view message,
          std::source_location location = std::source_location::current()) = 0;

    virtual void raw_log(std::string_view message) = 0;

    virtual void toggle_debug() = 0;

    virtual void toggle_print() = 0;

    virtual void toggle_alert() = 0;

    virtual void toggle_error() = 0;

    virtual void toggle_all() = 0;
  };

  template <typename... MsgType> std::string format(MsgType&&... msgs)
  {
    std::string result;
    (result += ... += std::forward<MsgType>(msgs));
    return result;
  }

} // namespace IO

#endif