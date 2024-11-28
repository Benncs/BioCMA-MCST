#ifndef __CLI_DYNLIB_LOADER_HPP__
#define __CLI_DYNLIB_LOADER_HPP__

#include <string_view>
#include <memory>

#include <dynlib/dynlib.hpp>

namespace UsafeUDF {
    [[nodiscard]] std::shared_ptr<DynamicLibrary> init_lib(std::string_view path);
}

#endif

