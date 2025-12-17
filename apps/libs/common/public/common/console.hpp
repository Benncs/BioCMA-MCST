#ifndef __LOGGER_CONSOLE_HPP__
#define __LOGGER_CONSOLE_HPP__

#include <common/logger.hpp>
#include <memory>
#include <optional>
#include <ostream>
#include <source_location>
#include <sstream>
#include <streambuf>

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

  // TODO move elsewhere

  enum RedirectionType
  {
    Buffer,
    File
  };

  class RedirectGuard;

  class RedirectHandle
  {
  public:
    RedirectHandle(std::shared_ptr<RedirectGuard> owner,
                   std::shared_ptr<bool> active_flag,
                   bool owns,
                   RedirectionType _type);

    ~RedirectHandle();

    RedirectHandle(const RedirectHandle&) = delete;
    RedirectHandle& operator=(const RedirectHandle&) = delete;

    RedirectHandle(RedirectHandle&&) = default;
    RedirectHandle& operator=(RedirectHandle&&) = default;

  private:
    std::shared_ptr<RedirectGuard> owner_ptr;
    std::shared_ptr<bool> active_flag;
    bool owns_guard = false;
    RedirectionType type;
  };

  class RedirectGuard : public std::enable_shared_from_this<RedirectGuard>
  {
  public:
    RedirectGuard();

    [[nodiscard("Handle must be held")]] RedirectHandle redirect();

    [[nodiscard("Handle must be held")]] RedirectHandle redirect_to_file();

    void restore();

    std::optional<std::string> getCapturedOutput() const;

  private:
    friend class RedirectHandle;

    bool has_been_redirected = false;
    std::streambuf* coutbuf = nullptr;
    std::stringstream buffer;
    int original_stdout_fd = -1;

    std::shared_ptr<bool> active_flag;
  };
} // namespace IO

#endif
