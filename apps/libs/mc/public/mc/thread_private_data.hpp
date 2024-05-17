#ifndef __MC_THREADS_PRIVATE_DATA_HPP__
#define __MC_THREADS_PRIVATE_DATA_HPP__

#include "mc/prng/prng.hpp"
#include <vector> 
#include <mc/particles/mcparticles.hpp>

namespace MC
{
  struct ThreadPrivateData
  {
    std::vector<MC::Particles> extra_process;
    std::vector<Particles *> in_dead_state;
    MC::PRNG rng;
  };
} // namespace MC


#endif //__MC__THREADS_PRIVATE_DATA_HPP__