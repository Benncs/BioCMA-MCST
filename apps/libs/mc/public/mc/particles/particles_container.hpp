#ifndef __PARTICLES_CONTAINER_HPP__
#define __PARTICLES_CONTAINER_HPP__

#include "Kokkos_ScatterView.hpp"
#include <cstddef>
#include <mc/particles/extra_process.hpp>
#include <mc/particles/particle_list.hpp>
#include <mc/particles/particle_model.hpp>

namespace MC
{
  /**
   * @brief Main owning object for Monte-Carlo particles.
   *
   * This class is responsible for managing and owning containers related to particles, providing
   * access to them for processing in a Monte-Carlo simulation. It uses a templated model to define
   * the behavior and properties of the particles.
   *
   * @tparam Model The particle model type used for defining the behavior and properties of the
   * particles.
   */
  template <ParticleModel Model> class ParticlesContainer
  {
  public:
    static constexpr std::size_t buffer_ratio = 1; // Buffer size = list.size()*buffer_ratio
    /**
     * @brief Alias for the model used by the container.
     */
    using UsedModel = Model;

    /**
     * @brief Default copy and move constructors and assignment operators.
     */
    ParticlesContainer() = default;
    ParticlesContainer(const ParticlesContainer&) = default;
    ParticlesContainer(ParticlesContainer&&) = default;
    ParticlesContainer& operator=(const ParticlesContainer&) = default;
    ParticlesContainer& operator=(ParticlesContainer&&) = default;

    /**
     * @brief Default destructor.
     */
    ~ParticlesContainer() = default;

    /**
     * @brief Returns a reference to the container holding the particles to be processed.
     *
     * @return Reference to the particle list in the compute space.
     */
    auto& get_compute() noexcept
    {
      return to_process;
    }

    /**
     * @brief Constructs a ParticlesContainer with a specified capacity.
     *
     * @param capacity The number of particles the container can initially hold.
     */
    explicit ParticlesContainer(size_t capacity) noexcept
        : to_process(capacity), host_process(), process_buffer(capacity * buffer_ratio, false)
    {
      // Initialization code can be placed here, if needed.
    }

    /**
     * @brief Gets the number of particles in the container.
     *
     * @return The number of particles currently in the to_process container.
     */
    [[nodiscard]] inline size_t n_particle() const noexcept
    {
      return to_process.size();
    }

    /**
     * @brief Migrates particles from compute space to host space and returns a reference to the
     * host container. A deep copy of main container is performed before returning.
     * @return Reference to the particle list in the host space.
     */
    auto& get_host()
    {
      migrate(to_process, host_process);
      return host_process;
    }

    /**
     * @brief Returns a reference to the extra results container.
     *
     * @return Reference to the results list in the compute space.
     */
    // auto& get_extra() noexcept
    // {
    //   return extra;
    // }

    template <class Archive> void save(Archive& ar) const
    {
      ar(to_process); //
    }

    template <class Archive> void load(Archive& ar)
    {
      ar(to_process); //
      process_buffer =
          ParticleList<ComputeSpace, Model>::with_capacity(to_process.size() * buffer_ratio);
    }

    MC::ParticleList<ComputeSpace, Model> process_buffer;

  private:
    ParticleList<ComputeSpace, Model> to_process; ///< Container for particles in the compute space.
    ParticleList<HostSpace, Model> host_process;  ///< Container for particles in the host space.
    // Results<ComputeSpace, Model> extra; ///< Container for additional results in the compute
    // space.
  };

} // namespace MC

#endif