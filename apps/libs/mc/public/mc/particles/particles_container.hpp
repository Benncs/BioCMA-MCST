#ifndef __PARTICLES_CONTAINER_HPP__
#define __PARTICLES_CONTAINER_HPP__


#include <Kokkos_Core.hpp>
#include <mc/traits.hpp>
#include <stdexcept>

namespace
{
  template <ModelType M> struct FillGapFunctor
  {

    FillGapFunctor(MC::ParticleStatus _status,
                   M::SelfParticle _model,
                   MC::ParticlePositions _position,
                   std::size_t _to_remove,
                   std::size_t _last_used_index)
        : status(std::move(_status)), model(std::move(_model)), position(std::move(_position)),
          offset("offset"), to_remove(_to_remove), last_used_index(_last_used_index)
    {

      Kokkos::deep_copy(offset, 0);
    }

    KOKKOS_INLINE_FUNCTION void operator()(const int i, std::size_t& update, const bool final)const 
    {
      const bool is_dead = (status(i) != MC::Status::Idle);

      std::size_t scan_index = update;
      update += is_dead ? 1 : 0;

      if (final && is_dead && scan_index < to_remove)
      {

        const auto i_to_remove = i;
        // Kokkos::printf("Removing index %d\r\n", i_to_remove);

        auto idx_to_move = last_used_index - Kokkos::atomic_fetch_add(&offset(), 1);
        while (status(idx_to_move) != MC::Status::Idle || idx_to_move == i_to_remove)
        {
          idx_to_move = last_used_index - Kokkos::atomic_fetch_add(&offset(), 1);
        }

        status(i_to_remove) = MC::Status::Idle;
        Kokkos::atomic_exchange(&position(i_to_remove), position(idx_to_move));

        for (std::size_t i_properties = 0; i_properties < M::n_var; ++i_properties)
        {
          model(i_to_remove, i_properties) = model(idx_to_move, i_properties);
        }
      }
    }

    MC::ParticleStatus status;
    M::SelfParticle model;
    MC::ParticlePositions position;
    Kokkos::View<std::size_t, ComputeSpace> offset;
    std::size_t to_remove;
    std::size_t last_used_index;
  };
}; // namespace

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
    /**
     * @brief Alias for the model used by the container.
     */
    using UsedModel = Model;

    explicit ParticlesContainer(std::size_t n_particle);
    ParticlesContainer();
    void init();

    Model::SelfParticle model;
    MC::ParticlePositions position;
    MC::ParticleStatus status;
    ParticleWeigths weights;
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
      return n_used_elements;
    };

    void get_contributions(std::size_t idx, const ContributionView& contributions) const;

    [[nodiscard]] KOKKOS_INLINE_FUNCTION bool handle_division(std::size_t idx1) const;

    void clean_dead(std::size_t to_remove);

    void merge_buffer();

    template <class Archive> void save(Archive& ar) const
    {
      (void)ar;
      // ar(to_process); //
      // process_buffer =
      //     ParticleList<ComputeSpace, Model>::with_capacity(to_process.size() * buffer_ratio);
    }

    [[nodiscard]] KOKKOS_INLINE_FUNCTION double get_weight(std::size_t idx) const;

  private:
    Model::SelfParticle buffer_model;
    ParticlePositions buffer_position;
    Kokkos::View<uint64_t, Kokkos::SharedSpace> buffer_index;
    double allocation_factor;
    std::size_t n_allocated_elements;
    uint64_t n_used_elements;

    static constexpr double default_allocation_factor = 2;
    void __allocate__(std::size_t new_size);
    void __allocate_buffer__();
    void __shrink__(std::size_t new_size, bool force);
  };

  template <ModelType Model>
  KOKKOS_INLINE_FUNCTION bool ParticlesContainer<Model>::handle_division(std::size_t idx1) const
  {
    PROFILE_SECTION("ParticlesContainer::handle_division")
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
    PROFILE_SECTION("ParticlesContainer::merge_buffer")
    const auto original_size = n_used_elements;
    const auto n_add_item = buffer_index();
    if (n_add_item == 0)
    {
      return;
    }
    __allocate__(original_size + n_add_item);
    // Merge position EZ
    Kokkos::parallel_for(
        "InsertNewPosition",
        Kokkos::RangePolicy<Kokkos::DefaultExecutionSpace>(0, n_add_item),
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
    PROFILE_SECTION("ParticlesContainer::__allocate__")
    if (new_size > 0)
    {
      if (new_size >= n_allocated_elements)
      {
        const auto new_allocated_size =
            static_cast<std::size_t>(std::ceil(static_cast<double>(new_size) * allocation_factor));
        n_allocated_elements = new_allocated_size;
        Kokkos::resize(position, n_allocated_elements);
        Kokkos::resize(model, n_allocated_elements);
        Kokkos::resize(status, n_allocated_elements);
        if constexpr (ConstWeightModelType<Model>)
        {
          Kokkos::resize(weights, n_allocated_elements);
        }
        else
        {
          Kokkos::resize(weights, 1);
        }
      }
    }
  }
  template <ModelType Model>
  void ParticlesContainer<Model>::__shrink__(std::size_t new_size, bool force)
  {
    PROFILE_SECTION("ParticlesContainer::__shrink__")
    if (new_size > 0 && (new_size > n_used_elements || force))
    {
      const auto new_allocated_size =
          static_cast<std::size_t>(std::ceil(static_cast<double>(new_size) * allocation_factor));
      n_allocated_elements = new_allocated_size;
      Kokkos::resize(position, n_allocated_elements);
      Kokkos::resize(model, n_allocated_elements);
      Kokkos::resize(status, n_allocated_elements);
      if constexpr (ConstWeightModelType<Model>)
      {
        Kokkos::resize(weights, n_allocated_elements);
      }
      else
      {
        Kokkos::resize(weights, 1);
      }
    }
  }

  template <ModelType Model> void ParticlesContainer<Model>::__allocate_buffer__()
  {
    PROFILE_SECTION("ParticlesContainer::__allocate_buffer__")
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

  template <ModelType M>
  ParticlesContainer<M>::ParticlesContainer(std::size_t n_particle)
      : model("particle_model", 0), position("particle_position", 0), // Dont allocate now
        weights("particle_weigth", 0),                                // Dont allocate now
        buffer_model("buffer_particle_model", 0),
        buffer_position("buffer_particle_position", 0), // Dont allocate now
        buffer_index("buffer_index"), allocation_factor(default_allocation_factor),
        n_allocated_elements(0), n_used_elements(n_particle)
  {
    __allocate__(n_particle);

    __allocate_buffer__();
  }

  template <ModelType M>
  ParticlesContainer<M>::ParticlesContainer()
      : model("particle_model"), position("particle_position"), weights("particle_weigth"),
        buffer_model("buffer_particle_model", 0),
        buffer_position("buffer_particle_position", 0), // Dont allocate now
        buffer_index("buffer_index"), allocation_factor(default_allocation_factor),
        n_allocated_elements(0), n_used_elements(0)

  {
  }

  template <ModelType M> void ParticlesContainer<M>::clean_dead(std::size_t to_remove)
  {
    PROFILE_SECTION("ParticlesContainer::remove_dead")
    if (to_remove == 0)
    {
      return;
    }

    if (to_remove == n_used_elements)
    {
      __shrink__(0, true); // TODO check this function if size=0
      n_used_elements = 0;
    }
    else if (to_remove > n_used_elements)
    {
      throw std::runtime_error(
          "clean_dead: Error in kernel cannot remove more element than existing");
    }
    else
    {

      const auto new_used_item = n_used_elements - to_remove;
      // Kokkos::printf("Sorting \r\n");
      // std::cout << "n_used_elements " << n_used_elements << std::endl;
      // std::cout << "To remove " << to_remove << std::endl;

      const auto last_used_index = n_used_elements - 1;
      Kokkos::parallel_scan("find_and_fill_gap",
                            Kokkos::RangePolicy<ComputeSpace>(0, n_used_elements),
                            FillGapFunctor<M>(status, model, position, to_remove, last_used_index));

      // Kokkos::View<std::size_t, ComputeSpace> offset("offset");
      // Kokkos::deep_copy(offset, 0);
      // Kokkos::parallel_scan(
      //     "find_and_fill_gap",
      //     Kokkos::RangePolicy<ComputeSpace>(0, n_used_elements),
      //     KOKKOS_CLASS_LAMBDA(const int i, std::size_t& update, const bool final) {
      //       const bool is_dead = (status(i) != Status::Idle);

      //       std::size_t scan_index = update;
      //       update += is_dead ? 1 : 0;

      //       if (final && is_dead && scan_index < to_remove)
      //       {

      //         const auto i_to_remove = i;
      //         // Kokkos::printf("Removing index %d\r\n", i_to_remove);

      //         auto idx_to_move = last_used_index - Kokkos::atomic_fetch_add(&offset(), 1);
      //         while (status(idx_to_move) != Status::Idle || idx_to_move == i_to_remove)
      //         {
      //           idx_to_move = last_used_index - Kokkos::atomic_fetch_add(&offset(), 1);
      //         }

      //         status(i_to_remove) = Status::Idle;
      //         Kokkos::atomic_exchange(&position(i_to_remove), position(idx_to_move));

      //         for (std::size_t i_properties = 0; i_properties < M::n_var; ++i_properties)
      //         {
      //           model(i_to_remove, i_properties) = model(idx_to_move, i_properties);
      //         }
      //       }
      //     });
      Kokkos::fence();
      KOKKOS_ASSERT(this->position.extent(0) == n_allocated_elements);
      KOKKOS_ASSERT(this->model.extent(0) == n_allocated_elements);
      KOKKOS_ASSERT(this->status.extent(0) == n_allocated_elements);
      n_used_elements = new_used_item;

      if (static_cast<double>(n_used_elements) / static_cast<double>(n_allocated_elements) <= 0.1)
      {
        // std::cout << "SHRINK" << std::endl;
        __shrink__(n_used_elements * 2, false);
      }
    };
  }

  template <ModelType M>
  void ParticlesContainer<M>::get_contributions(std::size_t idx,
                                                const ContributionView& contributions) const
  {
    static_assert(ConstWeightModelType<M>, "ModelType: Const apply_weight()");
    const double weight = get_weight(idx);
    M::contribution(idx, position(idx), weight, model, contributions);
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

} // namespace MC

#endif