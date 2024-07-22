#ifndef __MC_THREADS_PRIVATE_DATA_HPP__
#define __MC_THREADS_PRIVATE_DATA_HPP__

#include "common/execinfo.hpp"
#include "mc/prng/prng.hpp"
#include <cstdint>
#include <mc/particles/mcparticles.hpp>
#include <vector>

namespace MC
{
  struct alignas(ExecInfo::cache_line_size) ThreadPrivateData
  {
    std::vector<MC::Particles> extra_process;
    std::vector<std::int64_t> index_in_dead_state;
    MC::PRNG rng;
  };
} // namespace MC

#endif //__MC__THREADS_PRIVATE_DATA_HPP__