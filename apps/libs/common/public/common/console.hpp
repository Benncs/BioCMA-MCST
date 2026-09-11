#ifndef __LOGGER_CONSOLE_HPP__
#define __LOGGER_CONSOLE_HPP__

#include <common/logger.hpp>

#include <cstdint>
#include <ostream>
#include <source_location>

namespace IO
{

  /**
   * @brief Concrete class logging system through console.
   *
   * This class defines the interface for logging messages of various types
   */
  class Console final : public Logger
  {
  public:
    Console();
    ~Console() override = default;

    Console(const Console&) = delete;
    Console(Console&&) = delete;
    Console& operator=(const Console&) = delete;
    Console& operator=(Console&&) = delete;

    void debug(std::string_view message) noexcept final;

    void print(std::string_view prefix,
               std::string_view message) noexcept final;
    void alert(std::string_view prefix,
               std::string_view message) noexcept final;

    void error(std::string_view message,
               std::source_location location
               = std::source_location::current()) noexcept final;

    void raw_log(std::string_view message) noexcept final;

    void toggle_debug() noexcept final;
    void toggle_print() noexcept final;
    void toggle_alert() noexcept final;
    void toggle_error() noexcept final;
    void toggle_all() noexcept final;

    // std::ostream& get_output()
    // {
    //   return output;
    // };

  private:
    uint32_t flags;
    std::ostream& output;
    std::ostream& err_output;

    bool do_flush;
  };

} // namespace IO

#endif
