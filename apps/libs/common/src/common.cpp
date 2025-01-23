#include <biocma_cst_config.hpp>
#include <common/execinfo.hpp>
#include <string>

std::string ExecInfo::get_version()
{
  return std::string("v" + std::to_string(_BIOMC_VERSION_MAJOR) + "." +
                     std::to_string(_BIOMC_VERSION_MINOR) + "." +
                     std::to_string(_BIOMC_VERSION_DEV));
}
