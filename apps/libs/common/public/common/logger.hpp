#ifndef __COMMON_LOGGER_HPP__
#define __COMMON_LOGGER_HPP__

#include <memory>
#include <optional>
#include <sstream>
#include <streambuf>

namespace IO
{

  namespace AnsiCode
  {
    constexpr auto blue = "\033[34m";
    constexpr auto red = "\033[31m";
    constexpr auto reset = "\033[0m";

  } // namespace AnsiCode

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