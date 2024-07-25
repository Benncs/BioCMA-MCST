#ifndef __PARTICLES_CONTAINER_HPP__
#define __PARTICLES_CONTAINER_HPP__

#include "mc/particles/mcparticles.hpp"
#include <mc/particles/particles_list.hpp>
#include <mc/thread_private_data.hpp>

namespace MC
{

 
  class ParticlesContainer
  {

  public:
    explicit ParticlesContainer() = default;
    explicit ParticlesContainer(size_t capacity,
                                double weight) noexcept;
    explicit ParticlesContainer(const ParticlesContainer &other) = delete;
    ParticlesContainer(ParticlesContainer &&other) noexcept;
    ~ParticlesContainer() = default;

    ParticlesContainer &operator=(const ParticlesContainer &) = delete;
    ParticlesContainer &operator=(ParticlesContainer &&) = default;

    void merge(ThreadPrivateData& i_data);
    void init_extra(size_t n_extra);

    ParticlesList to_process;

    template <class Archive> void serialize(Archive &ar) { ar(to_process); }
  };

} // namespace MC

#endif