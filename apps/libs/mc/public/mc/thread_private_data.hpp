#ifndef __MC_THREADS_PRIVATE_DATA_HPP__
#define __MC_THREADS_PRIVATE_DATA_HPP__

#include <vector> 
#include <mc/particles/mcparticles.hpp>

namespace MC
{
  struct ThreadPrivateData
  {
    std::vector<MC::Particles> extra_process;
    std::vector<Particles *> in_dead_state;
  };
}


#endif //__MC__THREADS_PRIVATE_DATA_HPP__