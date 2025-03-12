#ifndef __PARTICLES_CONTAINER_HPP__
#define __PARTICLES_CONTAINER_HPP__

// #include "Kokkos_ScatterView.hpp"
// #include <cstddef>
// #include <mc/particles/extra_process.hpp>
// #include <mc/particles/particle_list.hpp>
// #include <mc/particles/particle_model.hpp>
#include <mc/traits.hpp>

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
  template <ModelType Model> class ParticlesContainer
  {
  public:
    static constexpr double buffer_ratio = 0.4; // Buffer size = ceil(list.size()*buffer_ratio)
    static constexpr std::size_t buffer_ratio = 1; // Buffer size = list.size()*buffer_ratio
    /**
     * @brief Alias for the model used by the container.
     */
    using UsedModel = Model;

    explicit ParticlesContainer(std::size_t n_particle);
    ParticlesContainer();
    void init();

    Model::SelfParticle model;
    MC::ParticlePositions position;

    /**
     * @brief Default copy and move constructors and assignment operators.
     */
    ParticlesContainer(const ParticlesContainer&) = default;
    ParticlesContainer(ParticlesContainer&&) = default;
    ParticlesContainer& operator=(const ParticlesContainer&) = default;
    ParticlesContainer& operator=(ParticlesContainer&&) = default;

    /**
     * @brief Default destructor.
     */
    ~ParticlesContainer() = default;

    /**
     * @brief Gets the number of particles in the container.
     *
     * @return The number of particles currently in the to_process container.
     */
    [[nodiscard]] KOKKOS_INLINE_FUNCTION std::size_t n_particles() const
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
      ParticleList<ComputeSpace, Model>::migrate(to_process, host_process);
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
      // ar(to_process); //
      // process_buffer =
      //     ParticleList<ComputeSpace, Model>::with_capacity(to_process.size() * buffer_ratio);
    }

  private:
    Model::SelfParticle buffer_model;
    ParticlePositions buffer_position;
    Kokkos::View<uint64_t, Kokkos::SharedSpace> buffer_index;
    void __allocate__(std::size_t new_size);
    void __allocate_buffer__();
    static constexpr double default_allocation_factor = 2;
    double allocation_factor;
    std::size_t n_allocated_elements;
    uint64_t n_used_elements;
  };

  template <ModelType Model>
  KOKKOS_INLINE_FUNCTION bool ParticlesContainer<Model>::handle_division(std::size_t idx1) const
  {
    if (Kokkos::atomic_load(&buffer_index()) < buffer_model.extent(0))
    {
      const auto idx2 = Kokkos::atomic_fetch_add(&buffer_index(), 1);
      Model::division(idx1, idx2, model, buffer_model);
      buffer_position(idx2) = position(idx1);
      return true;
    }
    return false;
  }

  template <ModelType Model> void ParticlesContainer<Model>::merge_buffer()
  {
    const auto original_size = n_used_elements;
    const auto n_add_item = buffer_index() + 1;
    __allocate__(original_size + n_add_item);
    // Merge position EZ
    Kokkos::parallel_for(
        "InsertNewPosition",
        Kokkos::RangePolicy<Kokkos::DefaultExecutionSpace>(0, n_add_item - 1),
        KOKKOS_CLASS_LAMBDA(const size_t i) { position(original_size + i) = buffer_position(i); });
    Kokkos::fence();
    // Merge model less EZ
    Kokkos::parallel_for(
        "InsertNewModel",
        Kokkos::RangePolicy<Kokkos::DefaultExecutionSpace>(0, n_add_item),
        KOKKOS_CLASS_LAMBDA(const size_t i) {
          for (std::size_t j = 0; j < Model::n_var; ++j)
          {
            model(original_size + i, j) = buffer_model(i, j);
          }
        });

    buffer_index() = 0;
    n_used_elements += n_add_item;
    __allocate_buffer__();
  }

  template <ModelType Model>
  void ParticlesContainer<Model>::__allocate__(const std::size_t new_size)
  {

    if (new_size > 0)
    {
      if (new_size >= n_allocated_elements)
      {
        const auto new_allocated_size =
            static_cast<std::size_t>(std::ceil(static_cast<double>(new_size) * allocation_factor));
        n_allocated_elements = new_allocated_size;
        Kokkos::resize(position, n_allocated_elements);
        Kokkos::resize(model, n_allocated_elements);
      }
    }
  }

  template <ModelType Model> void ParticlesContainer<Model>::__allocate_buffer__()
  {
    auto buffer_size = buffer_position.extent(0);
    if (static_cast<double>(buffer_size) / static_cast<double>(n_allocated_elements) < buffer_ratio)
    {
      buffer_size = std::ceil(static_cast<double>(n_allocated_elements) * buffer_ratio);
      // Realloc because not needed to keep buffer as it has been copied
      Kokkos::realloc(buffer_position, buffer_size);
      Kokkos::realloc(buffer_model, buffer_size);
      buffer_index() = 0;
    }
  }

} // namespace MC

#endif