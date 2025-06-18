#include "biocma_cst_config.hpp"
#include <common/logger.hpp>
#include <cstdio>
#include <iostream>
#include <optional>
#include <stdexcept>

constexpr bool f_redirect = true;

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

} // namespace IO