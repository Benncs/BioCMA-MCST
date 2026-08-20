#ifndef __CORE_RNG_SEED_HPP__
#define __CORE_RNG_SEED_HPP__

#include <common/execinfo.hpp>
#include <cstdint>
namespace Core
{
  uint64_t get_rng_seed([[maybe_unused]] const ExecInfo& info);
} // namespace Core
#endif
