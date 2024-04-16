#ifndef __PARTICLES_CONTAINER_HPP__
#define __PARTICLES_CONTAINER_HPP__

#include "mc/particles/mcparticles.hpp"
#include <mc/particles/particles_list.hpp>

namespace MC
{
  using STDVectorParticle = std::vector<MC::Particles>;

  struct TheadSafeData
  {
    STDVectorParticle extra_process;
    std::vector<Particles *> in_dead_state;
  };

  class ParticlesContainer
  {

  public:
    explicit ParticlesContainer() = default;
    explicit ParticlesContainer(size_t capacity, double weight) noexcept;
    explicit ParticlesContainer(const ParticlesContainer &other) = delete;
    explicit ParticlesContainer(ParticlesContainer &&other) noexcept;

    void merge(size_t i);
    void init_extra(size_t n_extra);

    ParticlesList to_process;
    std::vector<TheadSafeData> extras;
  };

} // namespace MC

#endif