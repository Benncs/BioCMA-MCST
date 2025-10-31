#ifndef __PARTICLES_CONTAINER_HPP__
#define __PARTICLES_CONTAINER_HPP__

#include "biocma_cst_config.hpp"
#include <Kokkos_Core.hpp>
#include <common/common.hpp>
#include <common/env_var.hpp>
#include <common/has_serialize.hpp>
#include <cstdint>
#include <mc/alias.hpp>
#include <mc/prng/prng.hpp>
#include <mc/traits.hpp>
namespace MC
{

  /**
   * @brief Main owning object for Monte-Carlo particles.
   *
   * This class is responsible for managing and owning containers related to
   * particles, providing access to them for processing in a Monte-Carlo
   * simulation. It uses a templated model to define the behavior and properties
   * of the particles.
   *
   * @tparam Model The particle model type used for defining the behavior and
   * properties of the particles.
   */
  template <ModelType Model> class ParticlesContainer
  {
  public:
    using UsedModel = Model;

    /**
     * @brief Alias for the model used by the container.
     */

    explicit ParticlesContainer(std::size_t n_particle);
    ParticlesContainer(); //=default;
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

    // NOLINTBEGIN(cppcoreguidelines-non-private-member-variables-in-classes)
    Model::SelfParticle model;
    MC::ParticlePositions position;
    MC::ParticleStatus status;
    ParticleWeigths weights;
    ParticleAges ages;
    // NOLINTEND(cppcoreguidelines-non-private-member-variables-in-classes)

    /**
     * @brief Get the contribution if particle at index idx
     * @tparam CviewType: view type fo contribution
     */
    template <typename CviewType>
    KOKKOS_INLINE_FUNCTION void
    get_contributions(const std::size_t idx,
                      const CviewType& contributions) const
    {

      static_assert(ConstWeightModelType<Model>,
                    "ModelType: Constapply_weight()");

      const double weight = get_weight(idx);
      const auto pos = position(idx);
      auto access = contributions.access();

      for (int i = begin; i < end; ++i)
      {
        const int rel = i - begin;
        const double c = weight * model(idx, i);
        access(rel, pos) += c;
      }
    }

    /**
     * @brief Attempts to spawn a new particle by performing division on the
     * specified particle.
     *
     * Executes container and buffer manipulations to spawn a new particle and
     * invokes `Model::division`. The operation may fail if conditions for
     * division are not met (e.g., insufficient space or invalid state).
     *
     * @param random_pool The random number generator pool for stochastic
     * operations.
     * @param idx1 The index of the particle to undergo division.
     * @return `true` if the spawn was successful, `false` otherwise.
     */
    [[nodiscard]] KOKKOS_INLINE_FUNCTION bool
    handle_division(const MC::KPRNG::pool_type& random_pool,
                    std::size_t idx1) const;

    /**
     * @brief Return the particle weight
     */
    [[nodiscard]] KOKKOS_INLINE_FUNCTION double
    get_weight(std::size_t idx) const;

    // HOST

    /**
     * @brief Clean all the non-idle particle even if number smaller than
     * threshold
     */
    void force_remove_dead();

    /**
     * @brief Insert particle buffer into the main container
     */
    void merge_buffer();

    /**
     * @brief Save data into ar for serialization
     */
    template <class Archive> void save(Archive& ar) const;

    /**
     * @brief Load data from ar for deserialization
     */
    template <class Archive> void load(Archive& ar);

    /**
     * @brief Get the allocation factor
     */
    [[nodiscard]] double get_allocation_factor() const noexcept;

    /**
     * @brief Returns the total number of elements the container can hold
     * without reallocating.
     */
    [[nodiscard]] std::size_t capacity() const noexcept;

    /**
     * @brief Returns the total number of non-idlde particle
     */
    [[nodiscard]] std::size_t get_inactive() const noexcept;

    /**
     * @brief Update the number of non-idle particle and remove if total number
     * exceeds threshold
     * May call : remove_inactive_particles
     */
    void update_and_remove_inactive(std::size_t out, std::size_t dead);

    /**
     * @brief Gets the number of particles in the container.
     *
     * @return The number of particles currently in the to_process
     * container.
     */
    [[nodiscard]] KOKKOS_INLINE_FUNCTION std::size_t n_particles() const;

    /**
     * @brief Removes inactive particles from the container.
     *
     * This function handles the removal of inactive particles based on the
     * specified count. Depending on the number of particles to remove, it may:
     *   - Shrink the container (i.e reallocate memory if necessary),
     *   - Do nothing if no removal is required,
     *   - Throw an error if `to_remove` is invalid (e.g., out of bounds or
     * incoherent),
     *   - Move an active particle to the position of the removed inactive
     * particle to reduce fragmentation.
     * @note: inactive_counter is updated at after remove
     *
     * @param to_remove The number of inactive particles to remove.
     * @throws std::runtime_error If `to_remove` is incoherent or exceeds valid
     * bounds.
     */
    void remove_inactive_particles(std::size_t to_remove);

// FIXME: used only in unit test
#ifndef NDEBUG
    [[maybe_unused]] [[nodiscard]] auto get_buffer_index() const;
#endif

  private:
    Model::SelfParticle buffer_model;
    ParticlePositions buffer_position;
    Kokkos::View<uint64_t, Kokkos::SharedSpace> buffer_index;
    std::size_t n_allocated_elements;
    uint64_t n_used_elements;
    std::size_t inactive_counter;
    void __allocate_buffer__();
    void _resize(std::size_t new_size, bool force = false);

    void load_tuning_constant();

    uint64_t minimum_dead_particle_removal{};
    double buffer_ratio{};
    double allocation_factor = {};
    double shink_ratio{};

    int begin;
    int end;
  };

} // namespace MC

/*
 * IMPLEMENTATION
 */
namespace MC
{
  namespace
  {

    using TeamPolicy = Kokkos::TeamPolicy<ComputeSpace>;
    using TeamMember = TeamPolicy::member_type;

    /**
     * @brief Functor for parallel particle removal and compaction using a
     * scan-like approach.
     *
     * This functor is designed to be executed in parallel (e.g., via Kokkos
     * parallel algorithms) to identify and remove non-idle particles
     * while compacting the container. It performs the following steps:
     *
     * 1. **Inactive Particle Detection**:
     *    - Checks if the particle at index `i` non-idle using
     *
     * 2. **Scan Phase**:
     *    - During the initial pass (`final = false`), it only counts inactive
     * particles.
     *    - The `scan_index` variable captures the prefix sum of inactive
     * particles up to index `i`.
     *
     * 3. **Compaction Phase** (final pass, `final = true`):
     *    - If the particle is  inactive and within the removal range
     * (`scan_index < to_remove`), it finds a replacement particle from the end
     * of the container (using atomic operations to avoid race conditions).
     *    - The replacement particle must be alive (idle) and not the same as
     * the particle being removed
     *
     * 4. **Particle Replacement**:
     *    - Copies the position, model properties, and age data from the
     * replacement particle to the position of the removed particle.
     *    - Marks the original position of the removed particle as idle.
     *
     * @note **Thread Safety**:
     *    - Uses atomic operations (`Kokkos::atomic_fetch_add`,
     * `Kokkos::atomic_exchange`) to ensure safe concurrent access to shared
     * variables (`update`, `offset`).
     *    - The `final` pass is executed exactly once per particle, so no atomic
     * operations are needed for the position and property copies.
     *
     *
     */
    template <ModelType M> struct CompactParticlesFunctor
    {

      CompactParticlesFunctor(MC::ParticleStatus _status,
                              M::SelfParticle _model,
                              MC::ParticlePositions _position,
                              MC::ParticleAges _ages,
                              std::size_t _to_remove,
                              std::size_t _last_used_index)
          : status(std::move(_status)), model(std::move(_model)),
            position(std::move(_position)), ages(std::move(_ages)),
            offset("offset"), to_remove(_to_remove),
            last_used_index(_last_used_index)
      {

        Kokkos::deep_copy(offset, 0);
      }

      KOKKOS_INLINE_FUNCTION void
      operator()(const int i, std::size_t& update, const bool final) const
      {

        /** @note **Performance Consideration**:
         *    - The loop to find a valid replacement particle (`idx_to_move`)
         * may introduce overhead if many particles are inactive or if the
         * container is highly fragmented.
         *    - TODO: Consider using hierarchical parallelism (e.g., thread
         * ranges) to optimize this step further.*/

        // Check if the particle is non-idle.
        const bool is_inactive = (status(i) != MC::Status::Idle);

        // Capture the current prefix sum of inactive particles.
        const std::size_t scan_index = update;

        // Increment the shared counter if the particle is inactive.
        update += is_inactive ? 1 : 0;

        // Final pass: Compact the container by replacing inactive particles.
        if (final && is_inactive && scan_index < to_remove)
        {

          const auto inactive_slot = i; // Index of the "gap" to fill

          // Atomically find a replacement particle from the end of the
          // container. The replacement must be alive (idle) and not the same as
          // the particle being removed.
          // Index of the particle that will "fill" the inactive slot
          auto replacement_index =
              last_used_index - Kokkos::atomic_fetch_add(&offset(), 1);
          while (status(replacement_index) != MC::Status::Idle ||
                 replacement_index == static_cast<std::size_t>(inactive_slot))
          {
            replacement_index =
                last_used_index - Kokkos::atomic_fetch_add(&offset(), 1);
          }

          // Mark the removed particle's position as idle.
          status(inactive_slot) = MC::Status::Idle;

          // Merge position EZ
          // No need atomic: final executed exactly once
          Kokkos::atomic_exchange(&position(inactive_slot),
                                  position(replacement_index));

          // TODO Use hierachical parallism here, thread range is likely to work
          for (std::size_t i_properties = 0; i_properties < M::n_var;
               ++i_properties)
          {
            model(inactive_slot, i_properties) =
                model(replacement_index, i_properties);
          }
          ages(inactive_slot, 0) = ages(replacement_index, 0);
          ages(inactive_slot, 1) = ages(replacement_index, 1);
        }
      }

      MC::ParticleStatus status;
      M::SelfParticle model;
      MC::ParticlePositions position;
      MC::ParticleAges ages;
      Kokkos::View<std::size_t, ComputeSpace> offset;
      std::size_t to_remove;
      std::size_t last_used_index;
    };

    /**
     * @brief Copies particle data from a buffer to the main container in
     * parallel.
     *
     * This functor is designed to be executed by a Kokkos team to:
     * 1. Copy model properties from buffer to main container with
     * back-insertion
     * 2. Copy the particle's position from `buffer_position` to the main
     * `position` container.
     * 3. Reset the particle's age values (hydraulic time) to 0 to avoid
     * additional buffer allocation.
     *
     *
     * @param team The Kokkos team member executing this functor.
     *             `team.league_rank()` determines the particle index `i`.
     */
    template <ModelType M> struct InsertFunctor
    {
      InsertFunctor(std::size_t _original_size,
                    M::SelfParticle _model,
                    MC::ParticlePositions _position,
                    MC::ParticleAges _ages,
                    M::SelfParticle _buffer_model,
                    MC::ParticlePositions _buffer_position)
          : original_size(_original_size), model(std::move(_model)),
            ages(std::move(_ages)), position(std::move(_position)),
            buffer_model(std::move(_buffer_model)),
            buffer_position(std::move(_buffer_position))
      {
      }

      // KOKKOS_INLINE_FUNCTION void operator()(const int i) const
      // {
      //   position(original_size + i) = buffer_position(i);
      //   for (std::size_t j = 0; j < M::n_var; ++j)
      //   {
      //     model(original_size + i, j) = buffer_model(i, j);
      //   }
      // }

      // TODO try this functor
      // TODO find a way to organise data to not copy non needed  data (like
      // contribs). Split model in two arrays?

      KOKKOS_INLINE_FUNCTION
      void operator()(const TeamMember& team) const
      {
        auto range = M::n_var;
        const int i = team.league_rank();

        Kokkos::parallel_for(
            Kokkos::TeamVectorRange(team, range),
            [&](const int& j)
            { model(original_size + i, j) = buffer_model(i, j); });
        position(original_size + i) = buffer_position(i);

        // Actually needs buffer to store mother's hydraulic time
        // But set new hydraulic time to 0 to not create new buffer a save
        // memory usage
        ages(original_size + i, 0) = 0;
        ages(original_size + i, 1) = 0;
      }

      std::size_t original_size;
      M::SelfParticle model;
      MC::ParticleAges ages;
      MC::ParticlePositions position;
      M::SelfParticle buffer_model;
      MC::ParticlePositions buffer_position;
    };

  }; // namespace

  template <ModelType Model>
  [[nodiscard]] KOKKOS_INLINE_FUNCTION std::size_t
  ParticlesContainer<Model>::n_particles() const
  {
    return n_used_elements;
  };

  template <ModelType Model>
  [[nodiscard]] std::size_t
  ParticlesContainer<Model>::get_inactive() const noexcept
  {
    return inactive_counter;
  }

  template <ModelType Model> void ParticlesContainer<Model>::force_remove_dead()
  {
    this->remove_inactive_particles(this->inactive_counter);
  }
#ifndef NDEBUG

  template <ModelType Model>
  [[maybe_unused]] [[nodiscard]] auto
  ParticlesContainer<Model>::get_buffer_index() const
  {
    return buffer_index();
  }
#endif

  template <ModelType Model>
  [[nodiscard]] double
  ParticlesContainer<Model>::get_allocation_factor() const noexcept
  {
    return allocation_factor;
  }

  template <ModelType Model>
  [[nodiscard]] std::size_t ParticlesContainer<Model>::capacity() const noexcept
  {
    return n_allocated_elements;
  }

  template <ModelType Model>
  template <class Archive>
  void ParticlesContainer<Model>::load(Archive& ar)
  {
    // Basically store everything that is usefull
    ar(n_allocated_elements, n_used_elements, allocation_factor);

    // Allocation is done in the deserialize
    deserialize_view(ar, weights);
    deserialize_view(ar, position);
    deserialize_view(ar, status);
    deserialize_view(ar, model);
    deserialize_view(ar, ages);
#ifndef NDEBUG
    Kokkos::printf("ParticlesContainer::load: Check if load_tuning_constant "
                   "works with different value");
#endif
    load_tuning_constant(); // TODO check if working
    __allocate_buffer__();  // Dont forget to allocate buffer
  }

  template <ModelType Model>
  template <class Archive>
  void ParticlesContainer<Model>::save(Archive& ar) const
  {

    // Basically store everything that is usefull
    ar(n_allocated_elements, n_used_elements, allocation_factor);
    serialize_view(ar, weights);
    serialize_view(ar, position);
    serialize_view(ar, status);
    serialize_view(ar, model);
    serialize_view(ar, ages);
  }

  template <ModelType Model>
  void
  ParticlesContainer<Model>::update_and_remove_inactive(const std::size_t out,
                                                        const std::size_t dead)
  {
    // Actually out particle and dead ones are just inactive
    inactive_counter += out;
    inactive_counter += dead;

    const auto _threshold =
        std::max(minimum_dead_particle_removal,
                 static_cast<uint64_t>(
                     static_cast<double>(n_used_elements) *
                     AutoGenerated::MC::dead_particle_ratio_threshold));

    // TODO: May change threshold as container is now cleaned before exporting,
    if (inactive_counter > _threshold)
    {
      this->remove_inactive_particles(inactive_counter);
    }
  }

  template <ModelType Model>
  KOKKOS_INLINE_FUNCTION bool ParticlesContainer<Model>::handle_division(
      const MC::KPRNG::pool_type& random_pool, std::size_t idx1) const
  {
    if (Kokkos::atomic_load(&buffer_index()) < buffer_model.extent(0))
    {
      const auto idx2 = Kokkos::atomic_fetch_add(&buffer_index(), 1);
      Model::division(random_pool, idx1, idx2, model, buffer_model);
      buffer_position(idx2) = position(idx1);
      ages(idx1, 1) = 0;
      return true;
    }
    return false;
  }

  template <ModelType Model> void ParticlesContainer<Model>::merge_buffer()
  {
    PROFILE_SECTION("ParticlesContainer::merge_buffer")
    const auto original_size = n_used_elements;
    const auto n_add_item = buffer_index();
    if (n_add_item == 0)
    {
      return;
    }
    _resize(original_size + n_add_item);
    Kokkos::parallel_for("insert_merge",
                         TeamPolicy(n_add_item, Kokkos::AUTO, Model::n_var),
                         InsertFunctor<Model>(original_size,
                                              model,
                                              position,
                                              ages,
                                              buffer_model,
                                              buffer_position));

    buffer_index() = 0;
    n_used_elements += n_add_item;
    __allocate_buffer__();
  }

  template <ModelType Model>
  void ParticlesContainer<Model>::_resize(std::size_t new_size, bool force)
  {
    PROFILE_SECTION("ParticlesContainer::_resize")

    // Ensure new_size is greater than zero
    if (new_size > 0)
    {
      // Determine if resizing is necessary based on the condition
      if (new_size > n_allocated_elements || force)
      {
        // Calculate the new allocated size
        const auto new_allocated_size = static_cast<std::size_t>(
            std::ceil(static_cast<double>(new_size) * allocation_factor));

        // Update the allocated size
        n_allocated_elements = new_allocated_size;

        // Perform the resizing on all relevant data containers
        Kokkos::resize(position, n_allocated_elements);
        Kokkos::resize(model,
                       n_allocated_elements,
                       Model::n_var); // use 2nd dim resize if dynamic
        Kokkos::resize(status, n_allocated_elements);
        Kokkos::resize(ages, n_allocated_elements);

        // Handle resizing for weights based on model type
        if constexpr (ConstWeightModelType<Model>)
        {
          Kokkos::resize(weights, 1); // Fixed size for ConstWeightModelType
        }
        else
        {
          Kokkos::resize(weights, n_allocated_elements);
        }
      }
    }
  }

  template <ModelType Model>
  void ParticlesContainer<Model>::__allocate_buffer__()
  {
    PROFILE_SECTION("ParticlesContainer::__allocate_buffer__")
    auto buffer_size = buffer_position.extent(0);
    if (static_cast<double>(buffer_size) /
            static_cast<double>(n_allocated_elements) <
        buffer_ratio)
    {
      buffer_size = static_cast<std::size_t>(
          std::ceil(static_cast<double>(n_allocated_elements) * buffer_ratio));

      // Realloc because not needed to keep buffer as it has been copied
      Kokkos::realloc(buffer_position, buffer_size);
      Kokkos::realloc(buffer_model,
                      buffer_size,
                      Model::n_var); // use 2nd dim resize if dynamic
      buffer_index() = 0;
    }
  }

  template <ModelType M>
  ParticlesContainer<M>::ParticlesContainer(std::size_t n_particle)
      : model(Kokkos::view_alloc(Kokkos::WithoutInitializing, "particle_model"),
              0),
        position(Kokkos::view_alloc(Kokkos::WithoutInitializing,
                                    "particle_position"),
                 0),
        status(
            Kokkos::view_alloc(Kokkos::WithoutInitializing, "particle_status"),
            0),
        weights(
            Kokkos::view_alloc(Kokkos::WithoutInitializing, "particle_weigth"),
            0),
        ages(Kokkos::view_alloc(Kokkos::WithoutInitializing, "particle_age"),
             0),
        buffer_model("buffer_particle_model", 0),
        buffer_position("buffer_particle_position", 0),
        buffer_index("buffer_index"), n_allocated_elements(0),
        n_used_elements(n_particle), inactive_counter(0)
  {
    load_tuning_constant();
    if (n_particle != 0)
    {
      _resize(n_particle);
      __allocate_buffer__();
    }

    auto bounds = M::get_bounds();
    begin = bounds.begin;
    end = bounds.end;
  }

  template <ModelType M>
  ParticlesContainer<M>::ParticlesContainer() : ParticlesContainer(0)
  {
  }

  template <ModelType M>
  void ParticlesContainer<M>::remove_inactive_particles(std::size_t to_remove)
  {

    PROFILE_SECTION("ParticlesContainer::remove_inactive_particles")
    if (to_remove == 0)
    {
      return;
    }

    if (to_remove == n_used_elements)
    {
      _resize(0, true);
      n_used_elements = 0;
      inactive_counter = 0;
    }
    else if (to_remove > n_used_elements)
    {
      throw std::runtime_error(
          "remove_inactive_particles: Error in kernel cannot remove more "
          "element than existing");
    }
    else
    {

      const auto new_used_item = n_used_elements - to_remove;

      const auto last_used_index = n_used_elements - 1;
      Kokkos::parallel_scan(
          "find_and_fill_gap",
          Kokkos::RangePolicy<ComputeSpace>(0, n_used_elements),
          CompactParticlesFunctor<M>(
              status, model, position, ages, to_remove, last_used_index));

      Kokkos::fence();
      KOKKOS_ASSERT(this->position.extent(0) == n_allocated_elements);
      KOKKOS_ASSERT(this->model.extent(0) == n_allocated_elements);
      KOKKOS_ASSERT(this->status.extent(0) == n_allocated_elements);
      n_used_elements = new_used_item;
      if (static_cast<double>(n_used_elements) /
              static_cast<double>(n_allocated_elements) <=
          shink_ratio)
      {
        _resize(n_used_elements * allocation_factor, false);
      }
      inactive_counter = inactive_counter - to_remove;
    };
  }

  template <ModelType M>
  [[nodiscard]] KOKKOS_INLINE_FUNCTION double
  ParticlesContainer<M>::get_weight(std::size_t idx) const
  {
    if constexpr (ConstWeightModelType<M>)
    {
      return weights(0);
    }
    else
    {
      return weights(idx);
    }
  }

  template <ModelType M> void ParticlesContainer<M>::load_tuning_constant()
  {
    auto get_valid_value_or =
        []<typename T>(
            std::string_view varname, T vdefault, T min_value, T max_value)
    {
      T tmp_value = Common::read_env_or(varname, vdefault);
      if (tmp_value > min_value && tmp_value <= max_value)
      {
        return tmp_value;
      }
      return vdefault;
    };

    minimum_dead_particle_removal = Common::read_env_or(
        "BIOMC_MC_MINIMUM_REMOVAL",
        AutoGenerated::MC::default_minimum_dead_particle_removal);

    buffer_ratio =
        get_valid_value_or("BIOMC_MC_BUFFER_RATIO",
                           AutoGenerated::MC::default_MC_buffer_ratio,
                           0.,
                           1.);

    constexpr double max_alloc = 5.;
    allocation_factor = get_valid_value_or(
        "BIOMC_MC_ALLOC_FACTOR",
        AutoGenerated::MC::default_particle_container_allocation_factor,
        0.,
        max_alloc);

    shink_ratio = get_valid_value_or("BIOMC_MC_SHRINK_RATIO",
                                     AutoGenerated::MC::default_shink_ratio,
                                     0.,
                                     1.);
  }
} // namespace MC

#endif
