#ifndef __COMMON_CONFIG_LOADER_HPP__
#define __COMMON_CONFIG_LOADER_HPP__

#include <memory>
#include <string_view>

namespace Common
{
  class ConfigFile
  {
  public:
  private:
    ConfigFile() = default;
  };

  class ConfigFactory
  {
  public:
    static std::shared_ptr<ConfigFile> get(std::string_view file_path);
  };

} // namespace Common

#endif