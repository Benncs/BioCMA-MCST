#ifndef __PARTICLES_CONTAINER_HPP__
#define __PARTICLES_CONTAINER_HPP__

#include <mc/particles/particles_list.hpp>

namespace MC
{
  class ParticlesContainer
  {
  public:
    ParticlesContainer() = default;
    ParticlesContainer(size_t capacity,double weight)
    {
      to_process = ParticlesList(capacity,weight);
    }

    ParticlesContainer(const ParticlesContainer &other) = delete;

    ParticlesContainer(ParticlesContainer &&other) noexcept
    {
      if (this != &other)
      {
        to_process = std::move(other.to_process);
      }
    }

    ParticlesList to_process;
    ParticlesList extra_process;
  };
} // namespace MC

#endif