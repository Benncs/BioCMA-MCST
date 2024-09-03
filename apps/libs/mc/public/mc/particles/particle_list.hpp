#ifndef __MC_PARTICLE_LIST_HPP__
#define __MC_PARTICLE_LIST_HPP__

#include <Kokkos_Core.hpp>
#include <common/kokkos_vector.hpp>
#include <mc/particles/mcparticles.hpp>
#include <mc/prng/prng.hpp>

#define __MC_MAX_PARTICLES_PER_PROCESS__ 100e6

namespace MC
{

  static constexpr size_t MAX_PARTICLE_BUFFER =
      __MC_MAX_PARTICLES_PER_PROCESS__;

  /**
   * @brief A container for managing a list of Monte-Carlo particles in a
   * specified memory space.
   *
   * This class extends the `KokkosVector` and provides additional functionality
   * for handling Monte-Carlo particles, such as spawning new particles and
   * initializing them with specific properties.
   *
   * @tparam MemorySpace The memory space in which the particles are stored
   * (e.g., device or host memory).
   * @tparam Model The particle model type used for defining the behavior and
   * properties of the particles.
   */
  template <typename MemorySpace, ParticleModel Model>
  class ParticleList final : public KokkosVector<Particle<Model>, MemorySpace>
  {
  public:
    /**
     * @brief Constructs a ParticleList with a given capacity.
     *
     * @param capacity The initial capacity of the particle list.
     * @param alloc Determines whether memory should be allocated upon
     * construction.
     */
    explicit ParticleList(size_t capacity, bool alloc = true)
        : KokkosVector<Particle<Model>, MemorySpace>(
              capacity, alloc, "particle_list")
    {
    }

    /**
     * @brief Creates a ParticleList with the specified capacity, without
     * allocating memory immediately.
     *
     * @param capacity The capacity of the particle list.
     * @return A ParticleList instance with the specified capacity.
     */
    auto static with_capacity(std::size_t capacity)
    {
      auto rhs = ParticleList(capacity, false);
      return rhs;
    }

    /**
     * @brief Default constructor for ParticleList.
     */
    ParticleList() = default;

    /**
     * @brief Spawns a new particle, initializes it, and adds it to the list.
     *
     * @return A pointer to the newly spawned particle, or nullptr if the
     * particle couldn't be created.
     */
    KOKKOS_FUNCTION Particle<Model> *spawn()
    {
      if (this->emplace(std::move(Particle<Model>())))
      {
        auto &p = this->back();
        p.init(rng_instance);
        return &p;
      }
      return nullptr;
    }

    /**
     * @brief Initializes the particles in the list with a given weight.
     *
     * This method sets initial properties for each particle, including
     * assigning a random container.
     *
     * @param weight The weight to be assigned to each particle.
     */
    void init(double weight)
    {
      auto local_data = this->data();
      KPRNG rng;
      Kokkos::parallel_for(
          "InitializeParticles",
          Kokkos::RangePolicy<Kokkos::DefaultExecutionSpace>(0, this->size()),
          KOKKOS_LAMBDA(const int i) {
            auto p = Particle<Model>(weight);
            p.properties.current_container = rng.uniform_u(0, 4);
            p.init(rng);
            local_data(i) = std::move(p);
          });

      Kokkos::fence();
    }

    /**
     * @brief Spawns and allocates a specified number of new particles with a
     * given weight.
     *
     * @param n The number of particles to spawn.
     * @param weight The weight to assign to the new particles.
     */
    void _spawn_alloc(size_t n, double weight)
    {
      ParticleList<MemorySpace, Model> newlist(n);
      newlist.init(weight);
      this->insert(newlist);
    }

    MC::KPRNG rng_instance; ///< Instance of the random number generator used
                            ///< for particle
  };

} // namespace MC

#endif