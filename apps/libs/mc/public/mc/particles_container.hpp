#ifndef __PARTICLES_CONTAINER_HPP__
#define __PARTICLES_CONTAINER_HPP__

#include "Kokkos_Core_fwd.hpp"
#include "common/common.hpp"
#include "mc/alias.hpp"
#include <Kokkos_Core.hpp>
#include <common/has_serialize.hpp>
#include <cstdint>
#include <mc/prng/prng.hpp>
#include <mc/traits.hpp>
#include <stdexcept>
#include <type_traits>

namespace
{

  using TeamPolicy = Kokkos::TeamPolicy<ComputeSpace>;
  using TeamMember = TeamPolicy::member_type;

  template <ModelType M> struct FillGapFunctor
  {

    FillGapFunctor(MC::ParticleStatus _status,
                   M::SelfParticle _model,
                   MC::ParticlePositions _position,
                   MC::ParticleAges _ages,
                   std::size_t _to_remove,
                   std::size_t _last_used_index)
        : status(std::move(_status)), model(std::move(_model)), position(std::move(_position)),
          ages(std::move(_ages)), offset("offset"), to_remove(_to_remove),
          last_used_index(_last_used_index)
    {

      Kokkos::deep_copy(offset, 0);
    }

    KOKKOS_INLINE_FUNCTION void operator()(const int i, std::size_t& update, const bool final) const
    {
      const bool is_dead = (status(i) != MC::Status::Idle);

      std::size_t scan_index = update;
      update += is_dead ? 1 : 0;

      if (final && is_dead && scan_index < to_remove)
      {

        const auto i_to_remove = i;
        auto idx_to_move = last_used_index - Kokkos::atomic_fetch_add(&offset(), 1);
        while (status(idx_to_move) != MC::Status::Idle ||
               idx_to_move == static_cast<std::size_t>(i_to_remove))
        {
          idx_to_move = last_used_index - Kokkos::atomic_fetch_add(&offset(), 1);
        }

        status(i_to_remove) = MC::Status::Idle;
        Kokkos::atomic_exchange(&position(i_to_remove), position(idx_to_move));
        // TODO Use hierachical parallism here, thread range is likely to work
        for (std::size_t i_properties = 0; i_properties < M::n_var; ++i_properties)
        {
          model(i_to_remove, i_properties) = model(idx_to_move, i_properties);
        }
        ages(i_to_remove, 0) = ages(idx_to_move, 0);
        ages(i_to_remove, 1) = ages(idx_to_move, 1);
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

  template <ModelType M> struct InsertFunctor
  {
    InsertFunctor(std::size_t _original_size,
                  M::SelfParticle _model,
                  MC::ParticlePositions _position,
                  M::SelfParticle _buffer_model,
                  MC::ParticlePositions _buffer_position)
        : original_size(_original_size), model(std::move(_model)), position(std::move(_position)),
          buffer_model(std::move(_buffer_model)), buffer_position(std::move(_buffer_position))
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
    // TODO find a way to organise data to not copy non needed  data (like contribs). Split model in
    // two arrays?

    KOKKOS_INLINE_FUNCTION
    void operator()(const TeamMember& team) const
    {
      auto range = M::n_var;
      const int i = team.league_rank();

      position(original_size + i) = buffer_position(i);
      Kokkos::parallel_for(Kokkos::TeamVectorRange(team, range),
                           [&](const int& j) { model(original_size + i, j) = buffer_model(i, j); });
    }

    std::size_t original_size;
    M::SelfParticle model;
    MC::ParticlePositions position;
    M::SelfParticle buffer_model;
    MC::ParticlePositions buffer_position;
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
    static constexpr double buffer_ratio =
        1; ///< Buffer size = ceil(list.size()*buffer_ratio), 0.6 is fine for high division rate or
           ///< stiff increase to 1
    /**
     * @brief Alias for the model used by the container.
     */
    using UsedModel = Model;

    explicit ParticlesContainer(std::size_t n_particle);
    ParticlesContainer(); //=default;

    Model::SelfParticle model;
    MC::ParticlePositions position;
    MC::ParticleStatus status;
    ParticleWeigths weights;
    ParticleAges ages;
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

    // TODO Add logic off internal_dead_counter header
    std::size_t counter = 0; //FIXME !!

    KOKKOS_INLINE_FUNCTION void get_contributions(std::size_t idx,
                                                  const ContributionView& contributions) const;

    [[nodiscard]] KOKKOS_INLINE_FUNCTION bool
    handle_division(const MC::KPRNG::pool_type& random_pool, std::size_t idx1) const;

    void clean_dead(std::size_t to_remove);

    void merge_buffer();

    [[maybe_unused]] [[nodiscard]] auto get_buffer_index() const
    {
      return buffer_index();
    }

    template <class Archive> void save(Archive& ar) const
    {
      ar(n_allocated_elements, n_used_elements, allocation_factor);
      serialize_view(ar, weights);
      serialize_view(ar, position);
      serialize_view(ar, status);
      serialize_view(ar, model);
      serialize_view(ar, ages);
    }

    template <class Archive> void load(Archive& ar)
    {
      ar(n_allocated_elements, n_used_elements, allocation_factor);
      deserialize_view(ar, weights);
      deserialize_view(ar, position);
      deserialize_view(ar, status);
      deserialize_view(ar, model);
      deserialize_view(ar, ages);

      __allocate_buffer__();
    }

    [[nodiscard]] KOKKOS_INLINE_FUNCTION double get_weight(std::size_t idx) const;

    [[nodiscard]] double get_allocation_factor() const
    {
      return allocation_factor;
    }

    [[nodiscard]] std::size_t capacity() const
    {
      return n_allocated_elements;
    }

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
  KOKKOS_INLINE_FUNCTION bool
  ParticlesContainer<Model>::handle_division(const MC::KPRNG::pool_type& random_pool,
                                             std::size_t idx1) const
  {
    if (Kokkos::atomic_load(&buffer_index()) < buffer_model.extent(0))
    {
      const auto idx2 = Kokkos::atomic_fetch_add(&buffer_index(), 1);
      Model::division(random_pool, idx1, idx2, model, buffer_model);
      buffer_position(idx2) = position(idx1);

      return true;
    }
    Kokkos::printf("%ld %ld\r\n", Kokkos::atomic_load(&buffer_index()), buffer_model.extent(0));
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

    // auto get_policy_insert = [=]()
    // {
    //   if constexpr (std::is_same_v<Kokkos::DefaultHostExecutionSpace,
    //                                Kokkos::DefaultExecutionSpace>)
    //   {
    //     return TeamPolicy(n_add_item / Kokkos::num_threads(), Kokkos::AUTO, Model::n_var);
    //   }
    //   else
    //   {
    //     return TeamPolicy(n_add_item, Kokkos::AUTO, Model::n_var);
    //   }
    // };

    Kokkos::parallel_for(
        "InsertMerge",
        TeamPolicy(n_add_item, Kokkos::AUTO, Model::n_var),
        InsertFunctor<Model>(original_size, model, position, buffer_model, buffer_position));

    buffer_index() = 0;
    n_used_elements += n_add_item;
    __allocate_buffer__();

  } 

  //TODO Merge reduce duplicate __allocate__ and __shrink__
  template <ModelType Model>
  void ParticlesContainer<Model>::__allocate__(const std::size_t new_size)
  {
    PROFILE_SECTION("ParticlesContainer::__allocate__")
    if (new_size > 0)
    {
      if (new_size >= this->n_allocated_elements)
      {
        const auto new_allocated_size =
            static_cast<std::size_t>(std::ceil(static_cast<double>(new_size) * allocation_factor));
        this->n_allocated_elements = new_allocated_size;
        Kokkos::resize(position, n_allocated_elements);
        Kokkos::resize(model, n_allocated_elements, Model::n_var); // use 2nd dim resize if dynamic
        Kokkos::resize(status, n_allocated_elements);
        Kokkos::resize(ages, n_allocated_elements);
        if constexpr (ConstWeightModelType<Model>)
        {
          Kokkos::resize(weights, 1);
        }
        else
        {
          Kokkos::resize(weights, n_allocated_elements);
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
      Kokkos::resize(model, n_allocated_elements, Model::n_var); // use 2nd dim resize if dynamic
      Kokkos::resize(status, n_allocated_elements);
      Kokkos::resize(ages, n_allocated_elements);
      if constexpr (ConstWeightModelType<Model>)
      {
        Kokkos::resize(weights, 1);
      }
      else
      {
        Kokkos::resize(weights, n_allocated_elements);
      }
    }
  }

  template <ModelType Model> void ParticlesContainer<Model>::__allocate_buffer__()
  {
    PROFILE_SECTION("ParticlesContainer::__allocate_buffer__")
    auto buffer_size = buffer_position.extent(0);
    if (static_cast<double>(buffer_size) / static_cast<double>(n_allocated_elements) < buffer_ratio)
    {
      buffer_size = static_cast<std::size_t>(
          std::ceil(static_cast<double>(n_allocated_elements) * buffer_ratio));

      // Realloc because not needed to keep buffer as it has been copied
      Kokkos::realloc(buffer_position, buffer_size);
      Kokkos::realloc(buffer_model, buffer_size, Model::n_var); // use 2nd dim resize if dynamic
      buffer_index() = 0;
    }
  }

  template <ModelType M>
  ParticlesContainer<M>::ParticlesContainer(std::size_t n_particle)
      : model(Kokkos::view_alloc(Kokkos::WithoutInitializing, "particle_model"), 0),
        position(Kokkos::view_alloc(Kokkos::WithoutInitializing, "particle_position"), 0),
        status(Kokkos::view_alloc(Kokkos::WithoutInitializing, "particle_status"), 0),
        weights(Kokkos::view_alloc(Kokkos::WithoutInitializing, "particle_weigth"), 0),
        ages(Kokkos::view_alloc(Kokkos::WithoutInitializing, "particle_age"), 0),
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
      : model(Kokkos::view_alloc(Kokkos::WithoutInitializing, "particle_model"), 0),
        position(Kokkos::view_alloc(Kokkos::WithoutInitializing, "particle_position"), 0),
        status(Kokkos::view_alloc(Kokkos::WithoutInitializing, "particle_status"), 0),
        weights(Kokkos::view_alloc(Kokkos::WithoutInitializing, "particle_weigth"), 0),
        ages(Kokkos::view_alloc(Kokkos::WithoutInitializing, "particle_age"), 0),
        buffer_model(Kokkos::view_alloc(Kokkos::WithoutInitializing, "buffer_particle_model"), 0),
        buffer_position(Kokkos::view_alloc(Kokkos::WithoutInitializing, "buffer_particle_model")),
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
      __shrink__(0, true);
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

      const auto last_used_index = n_used_elements - 1;
      Kokkos::parallel_scan(
          "find_and_fill_gap",
          Kokkos::RangePolicy<ComputeSpace>(0, n_used_elements),
          FillGapFunctor<M>(status, model, position, ages, to_remove, last_used_index));

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
  KOKKOS_INLINE_FUNCTION void
  ParticlesContainer<M>::get_contributions(std::size_t idx,
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