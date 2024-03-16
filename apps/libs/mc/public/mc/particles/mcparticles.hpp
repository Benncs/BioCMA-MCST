#ifndef __MC_PARTICLESHPP__
#define __MC_PARTICLESHPP__

#include <any>
#include <cstddef>
namespace MC
{
  class Particles
  {
  public:
    Particles() = default;
    size_t current_container;
    size_t current_domain;
    size_t random_seed;
    size_t id;
    std::any data;
  };
} // namespace MC

#endif