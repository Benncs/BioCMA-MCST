#ifndef __PARTICLES_CONTAINER_HPP__
#define __PARTICLES_CONTAINER_HPP__

#include <mc/particles/particles_list.hpp>

namespace MC
{
  class ParticlesContainer
  {
  public:
    ParticlesContainer() = default;
    ParticlesContainer(size_t capacity)
    {
      to_process = ParticlesList(capacity);
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
  };
} // namespace MC

#endif