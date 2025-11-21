#ifndef __UDF_HANDLE_HPP__
#define __UDF_HANDLE_HPP__

#include <common/results.hpp>
#include <memory>
#ifdef USE_UDF
#  include <udf_includes.hpp>
#endif

namespace Unsafe
{
#ifndef USE_UDF
  struct DynamicLibrary;
#endif
  Result<std::shared_ptr<DynamicLibrary>, std::string>
  load_udf(std::string_view model_name);
} // namespace Unsafe
#endif
