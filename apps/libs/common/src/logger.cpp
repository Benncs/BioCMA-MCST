#include <biocma_cst_config.hpp>
#include <common/logger.hpp>
#include <cstdio>
#include <iostream>
#include <optional>
#include <stdexcept>
#include <memory>
#include <optional>
#include <ostream>
#include <source_location>
#include <sstream>
#include <streambuf>
#include <common/console.hpp>

#include <unistd.h> //dup2

constexpr bool f_redirect = true;

namespace
{
  // NOLINTBEGIN(hicpp-signed-bitwise)
  enum Flags
  {
    None = 0,
    Debug = 1 << 0,
    Print = 1 << 1,
    Alert = 1 << 2,
    Error = 1 << 3,
    All = Debug | Print | Alert | Error
  };
  // NOLINTEND(hicpp-signed-bitwise)
} // namespace

namespace IO
{

  RedirectHandle::RedirectHandle(std::shared_ptr<RedirectGuard> owner,
                                 std::shared_ptr<bool> active_flag,
                                 bool owns,
                                 RedirectionType _type)
      : owner_ptr(std::move(owner)), active_flag(std::move(active_flag)), owns_guard(owns),
        type(_type)
  {
  }

  RedirectGuard::RedirectGuard() : active_flag(std::make_shared<bool>(false))
  {
  }

  RedirectHandle RedirectGuard::redirect()
  {

    if (*active_flag)
    {
      return {shared_from_this(), active_flag, false, RedirectionType::Buffer};
    }

    *active_flag = true;
    has_been_redirected = true;

    coutbuf = std::cout.rdbuf();

    std::cout.rdbuf(buffer.rdbuf());

    original_stdout_fd = dup(fileno(stdout));
    if (original_stdout_fd == -1)
    {
      throw std::runtime_error("dup failed");
    }

    return {shared_from_this(), active_flag, true, RedirectionType::Buffer};
  }

  RedirectHandle RedirectGuard::redirect_to_file()
  {
    if (*active_flag)
    {
      return {shared_from_this(), active_flag, false, RedirectionType::File};
    }

    *active_flag = true;
    has_been_redirected = false;

    coutbuf = std::cout.rdbuf();

    std::cout.rdbuf(buffer.rdbuf());

    original_stdout_fd = dup(fileno(stdout));
    if (original_stdout_fd == -1)
    {
      throw std::runtime_error("dup failed");
    }

    if (!freopen("log", "w", stdout))
    {
      throw std::runtime_error("freopen failed");
    }
    if (coutbuf != nullptr)
    {
      std::cout.rdbuf(coutbuf);
    }

    return {shared_from_this(), active_flag, true, RedirectionType::File};
  }

  void RedirectGuard::restore()
  {
    if (!*active_flag)
    {
      return;
    }

    if (coutbuf != nullptr)
    {
      std::cout.rdbuf(coutbuf);
      coutbuf = nullptr;
    }

    if (original_stdout_fd != -1)
    {
      fflush(stdout);
      if (dup2(original_stdout_fd, fileno(stdout)) == -1)
      {
        perror("dup2 restore");
      }
      close(original_stdout_fd);
      original_stdout_fd = -1;
    }

    *active_flag = false;
  }

  std::optional<std::string> RedirectGuard::getCapturedOutput() const
  {
    if (has_been_redirected)
    {
      return buffer.str();
    }
    return std::nullopt;
  }

  // RAII Restore logic
  RedirectHandle::~RedirectHandle()
  {
    if (owns_guard && active_flag && *active_flag)
    {
      owner_ptr->restore();
    }
  }

  Console::Console() : flags(Flags::None), output(std::cout),err_output(std::cerr)
  {
      std::ios::sync_with_stdio(false);
      std::locale::global(std::locale("en_US.utf8"));
  }

  void Console::debug([[maybe_unused]] std::string_view message)
  {
#ifndef NDEBUG
    if ((flags & Flags::Debug) != 0U)
    {
      output << "[DEBUG]: " << message << "\r\n";
    }
#endif
  }

  void Console::print(std::string_view prefix, std::string_view message)
  {
    if ((flags & Flags::Print) != 0U)
    {
      output <<AnsiCode::blue<<"["<<prefix<<"]: "<<AnsiCode::reset << message << "\r\n";
    }
  }

  void Console::alert(std::string_view prefix, std::string_view message)
  {
    if ((flags & Flags::Alert) != 0U)
    {
      output << AnsiCode::red<<"[Warning: "<<prefix<<"]" << ": "<<AnsiCode::reset  << message << "\r\n";
    }
  }

  void Console::error(std::string_view message,std::source_location location )
  {
    if ((flags & Flags::Error) != 0U)
    {
      err_output << AnsiCode::red<<"[Error] From: "<<location.function_name() <<AnsiCode::reset  <<": "<< message << "\r\n";
    }
  }

   void Console::raw_log(std::string_view message )
  {
    if ((flags & Flags::Error) != 0U)
    {
      err_output << message;
    }
  }

  void Console::toggle_debug()
  {

    flags ^= Flags::Debug;
  }

  void Console::toggle_all()
  {
    flags ^= Flags::All;
  }


  void Console::toggle_print()
  {
    flags ^= Flags::Print;
  }

  void Console::toggle_alert()
  {
    flags ^= Flags::Alert;
  }

  void Console::toggle_error()
  {
    flags ^= Flags::Error;
  }

} // namespace IO