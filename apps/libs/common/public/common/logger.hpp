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

  /**
   * @brief Base class for a generic logging system.
   *
   * This class defines the interface for logging messages of various types
   */
  class Logger
  {
  public:
    Logger() = default;
    Logger(const Logger&) = default;
    Logger(Logger&&) = default;
    Logger& operator=(const Logger&) = default;
    Logger& operator=(Logger&&) = default;
    virtual ~Logger() = default;

    /**
     * @brief Log a debug message.
     *
     * Logs a debug message with no prefix.
     *
     * @param message The message to log.
     */
    virtual void debug(std::string_view message) = 0;

    /**
     * @brief Log a regular print message with a prefix.
     *
     * Logs a general message with the specified prefix. Typically used for
     * general informational messages.
     *
     * @param prefix A string prefix to prepend to the message.
     * @param message The message to log.
     */
    virtual void print(std::string_view prefix, std::string_view message) = 0;

    /**
     * @brief Log an alert message with a prefix.
     *
     * Logs an alert message that typically requires immediate attention.
     *
     * @param prefix A string prefix to prepend to the message.
     * @param message The message to log.
     */
    virtual void alert(std::string_view prefix, std::string_view message) = 0;

    /**
     * @brief Log an error message with optional location details.
     *
     * Logs an error message, optionally including the location (file, line,
     * function) where the error occurred.
     *
     * @param message The error message to log.
     * @param location The source location of the error (default: current
     * location).
     */
    virtual void
    error(std::string_view message,
          std::source_location location = std::source_location::current()) = 0;

    /**
     * @brief Log a raw message with no additional formatting.
     *
     * Logs a raw message without any prefix or other enhancements.
     *
     * @param message The raw message to log.
     */
    virtual void raw_log(std::string_view message) = 0;

    /**
     * @brief Toggle the debug logging behavior on or off.
     *
     * Enables or disables the logging of debug messages.
     */
    virtual void toggle_debug() noexcept = 0;

    /**
     * @brief Toggle the print logging behavior on or off.
     *
     * Enables or disables the logging of print messages.
     */
    virtual void toggle_print() noexcept = 0;

    /**
     * @brief Toggle the alert logging behavior on or off.
     *
     * Enables or disables the logging of alert messages.
     */
    virtual void toggle_alert() noexcept = 0;

    /**
     * @brief Toggle the error logging behavior on or off.
     *
     * Enables or disables the logging of error messages.
     */
    virtual void toggle_error() noexcept = 0;

    /**
     * @brief Toggle all logging types on or off.
     *
     * Enables or disables all types of logging at once (debug, print, alert,
     * error).
     */
    virtual void toggle_all() noexcept = 0;
  };

  template <typename... MsgType> std::string format(MsgType&&... msgs)
  {
    std::string result;
    (result += ... += std::forward<MsgType>(msgs));
    return result;
  }

} // namespace IO

#endif