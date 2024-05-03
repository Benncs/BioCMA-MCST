#ifndef __PARTICLES_CONTAINER_HPP__
#define __PARTICLES_CONTAINER_HPP__

#include "mc/particles/mcparticles.hpp"
#include <mc/particles/particles_list.hpp>

namespace MC
{

  struct TheadSafeData
  {
    std::vector<MC::Particles> extra_process;
    std::vector<Particles *> in_dead_state;
  };

  class ParticlesContainer
  {

  public:
    explicit ParticlesContainer() = default;
    explicit ParticlesContainer(size_t capacity,
                                double weight,
                                size_t n_extra) noexcept;
    explicit ParticlesContainer(const ParticlesContainer &other) = delete;
    ParticlesContainer(ParticlesContainer &&other) noexcept;
    ~ParticlesContainer() = default;

    ParticlesContainer &operator=(const ParticlesContainer &) = delete;
    ParticlesContainer &operator=(ParticlesContainer &&) = default;

    void merge(size_t i);
    void init_extra(size_t n_extra);

    ParticlesList to_process;
    std::vector<TheadSafeData> extras; //TODO move to unit 
  };

} // namespace MC

#endif