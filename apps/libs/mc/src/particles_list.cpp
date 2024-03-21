#include "mc/particles/mcparticles.hpp"
#include <mc/particles/particles_list.hpp>

namespace MC
{
  // ParticlesList::ParticlesList(const ParticlesList &other) : data(other.data)
  // {
  // }

  ParticlesList::ParticlesList(size_t capacity,double weight)
  {
    this->data = std::vector<MC::Particles>(capacity,Particles(weight));
  }

  ParticlesList::ParticlesList(ParticlesList &&other) noexcept
      : data(std::move(other.data))
  {
  }

  // ParticlesList &ParticlesList::operator=(const ParticlesList &other)
  // {
  //   if (this != &other)
  //   {
  //     data = other.data;
  //   }
  //   return *this;
  // }
  ParticlesList &ParticlesList::operator=(ParticlesList &&other) noexcept
  {
    if (this != &other)
    {
      data = std::move(other.data);
    }
    return *this;
  }
} // namespace MC