#ifndef __MC_PARTICLE_LIST_HPP__
#define __MC_PARTICLE_LIST_HPP__

#include "common/common.hpp"
#include "mc/particles/data_holder.hpp"
#include <Kokkos_Assert.hpp>
#include <Kokkos_Core.hpp>
#include <Kokkos_Macros.hpp>
#include <biocma_cst_config.hpp>
#include <common/kokkos_vector.hpp>
#include <cstddef>
#include <cstdint>
#include <decl/Kokkos_Declare_OPENMP.hpp>
#include <mc/particles/mcparticles.hpp>
#include <mc/prng/prng.hpp>
#include <stdexcept>
#include <string>

#include <Kokkos_StdAlgorithms.hpp>



WARN_EXPERIMENTAL

template <typename Space, typename ViewType, class UnaryPredicateType>
void remove_if(uint64_t i_max_used,
               uint64_t remove_count,
               ViewType view,
               UnaryPredicateType pred)
{
  
  const size_t view_size = view.size();
  if (remove_count >= view_size)
  {
    return;
  }

  Space space;

  const uint64_t to_keep = i_max_used - remove_count;

  // Initialisation, make tmp view as keeping storage
  ViewType tmp_view(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, space, "remove_tmp_view"),
      view_size);

  uint64_t scan_count = 0;
  Kokkos::parallel_scan(
      "remove_if_scan",
      Kokkos::RangePolicy<Space>(space, 0, i_max_used),
      KOKKOS_LAMBDA(const auto &i, auto &update, const bool is_final) {
        const auto &val = view[i];
        if (!pred(val))
        {
          if (is_final)
          {
            tmp_view[update] = std::move(val);
          }
          update++;
        }
      },
      scan_count);

  Kokkos::fence("scan:remove_if");
  KOKKOS_ASSERT(scan_count == to_keep);
  Kokkos::parallel_for(
      "remove_if_parfor",
      Kokkos::RangePolicy<Space>(space, 0, i_max_used),
      KOKKOS_LAMBDA(const auto i) { view[i] = std::move(tmp_view[i]); });
  Kokkos::fence("parfor:remove_if");
  // view = tmp_view;
}

namespace MC
{

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
      if (capacity > AutoGenerated::MC_MAX_PARTICLE_BUFFER)
      {
        std::string err =
            "Cannot allocated more than " +
            std::to_string(AutoGenerated::MC_MAX_PARTICLE_BUFFER) +
            " particles";
        throw std::runtime_error(err);
      }
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

    void remove_dead(uint64_t count_dead)
    {
      PROFILE_SECTION("remove_dead")

      auto pred = KOKKOS_LAMBDA(const auto &p)
      {
        return p.properties.status == MC::CellStatus::DEAD ||
               p.properties.status == MC::CellStatus::OUT;
      };
      // auto beg = Kokkos::Experimental::begin(this->_owned_data);
      // auto end = beg + this->size();
      uint64_t new_used_item = this->size() - count_dead;
      remove_if<Kokkos::DefaultExecutionSpace>(
          this->size(), count_dead, this->_owned_data, pred);

      // auto res = Kokkos::Experimental::remove_if(
      //     "dead_removing", Kokkos::DefaultExecutionSpace(), beg, end, pred);

      // const auto stepsA = Kokkos::Experimental::distance(
      //     Kokkos::Experimental::begin(this->_owned_data), res);

      KOKKOS_ASSERT(this->_owned_data.extent(0) == this->capacity());
      this->set_n_used_elements(new_used_item);
      KOKKOS_ASSERT(this->size() <= this->capacity());
    }

    MC::KPRNG rng_instance; ///< Instance of the random number generator used
                            ///< for particle
  };

} // namespace MC

#endif