#ifndef __BIOMC_COMMON_ENV_VAR_HPP__
#define __BIOMC_COMMON_ENV_VAR_HPP__
#include <cstdio>
#include <cstdlib>
#include <optional>
#include <ostream>
#include <sstream>
#include <string_view>
namespace Common
{

  template <typename T> std::optional<T> read_env(std::string_view varname)
  {
    const char* env_var = std::getenv(varname.data());
    std::optional<T> ret = std::nullopt;
    if (env_var != nullptr)
    {
      std::istringstream stream(env_var);
      T value;
      stream >> value;
      if (stream)
      {
        std::printf(
            "[Config] Env value load: %s %s\r\n", varname.data(), env_var);
        std::fflush(stdout);

        ret = value;
      }
    }
    return ret;
  }

  /**
   * @brief Wrapper arround get_env to get envariable with fallback to default
   * @note: UB if T is not trivial type
   */
  template <typename T> T read_env_or(std::string_view varname, T vdefault)
  {
    return read_env<T>(varname).value_or(vdefault);
  }

  /**
   * @brief Wrapper arround set_env
   * @note: UB if T is not trivial type
   */
  template <typename T> bool set_local_env(std::string_view varname, T val)
  {
#ifdef _WIN32
    // Windows-specific handling (not implemented)
    std::cerr << "Not implemented for Win32" << std::endl;
    return false;
#else
    std::ostringstream oss;
    oss << val;
    std::string str_value = oss.str();

    return setenv(varname.data(), str_value.c_str(), 1) == 0;
#endif
  }

} // namespace Common

#endif
