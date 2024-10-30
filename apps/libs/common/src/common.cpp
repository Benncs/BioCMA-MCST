#include "common/execinfo.hpp"
#include <common/common.hpp>

std::string ExecInfo::get_version()
{
  return std::string("v" + std::to_string(_BIOMC_VERSION_MAJOR) + "." + std::to_string(_BIOMC_VERSION_MINOR) + "." +
                     std::to_string(_BIOMC_VERSION_DEV));
}