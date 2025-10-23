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
  class Console final : public Logger
  {
  public:
    Console();
    ~Console() override = default;

    Console(const Console&) = delete;
    Console(Console&&) = delete;
    Console& operator=(const Console&) = delete;
    Console& operator=(Console&&) = delete;

    void debug(std::string_view message) final;
    void print(std::string_view prefix, std::string_view message) final;
    void alert(std::string_view prefix, std::string_view message) final;
    void error(
        std::string_view message,
        std::source_location location = std::source_location::current()) final;

    void raw_log(std::string_view message) final;

    void toggle_debug() final;
    void toggle_print() final;
    void toggle_alert() final;
    void toggle_error() final;
    void toggle_all() final;

    std::ostream& get_output()
    {
      return output;
    };

  private:
    uint32_t flags;
    std::ostream& output;
    std::ostream& err_output;
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