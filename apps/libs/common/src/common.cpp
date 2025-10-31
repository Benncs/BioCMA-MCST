#include <biocma_cst_config.hpp>
#include <common/env_var.hpp>
#include <common/execinfo.hpp>
#include <string>
std::string ExecInfo::get_version()
{
  return std::string("v" + std::to_string(_BIOMC_VERSION_MAJOR) + "." +
                     std::to_string(_BIOMC_VERSION_MINOR) + "." +
                     std::to_string(_BIOMC_VERSION_DEV));
}
namespace Common
{
  template <> bool read_env_or(std::string_view varname, bool vdefault)
  {
    const char* env_var = std::getenv(varname.data());
    if (env_var != nullptr)
    {
      if (std::string(env_var) == "true" || std::string(env_var) == "1")
      {
        return true;
      }
      else if (std::string(env_var) == "false" || std::string(env_var) == "0")
      {
        return false;
      }
    }
    return vdefault;
  }

  template <>
  bool set_local_env<std::string>(std::string_view varname, std::string val)
  {
#ifdef _WIN32
    std::cerr << "Not implemented for Win32" << std::endl;
    return false;
#else
    return setenv(varname.data(), val.c_str(), 1) == 0;
#endif
  }

  template <> bool set_local_env<bool>(std::string_view varname, bool val)
  {
#ifdef _WIN32
    std::cerr << "Not implemented for Win32" << std::endl;
    return false;
#else
    return setenv(varname.data(), val ? "true" : "false", 1) == 0;
#endif
  }

} // namespace Common
