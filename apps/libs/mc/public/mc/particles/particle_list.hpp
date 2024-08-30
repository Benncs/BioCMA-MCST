#ifndef __MC_PARTICLE_LIST_HPP__
#define __MC_PARTICLE_LIST_HPP__

#include <Kokkos_Core_fwd.hpp>
#include <Kokkos_DynamicView.hpp>
#include <common/kokkos_vector.hpp>
#include <mc/particles/mcparticles.hpp>
#include <mc/prng/prng.hpp>

#define __MC_MAX_PARTICLES_PER_PROCESS__ 100e6

namespace MC
{

  static constexpr size_t MAX_PARTICLE_BUFFER =
      __MC_MAX_PARTICLES_PER_PROCESS__;

  template <typename MemorySpace, ParticleModel Model>
  class ParticleList final : public KokkosVector<Particle<Model>, MemorySpace>
  {
  public:
    explicit ParticleList(size_t capacity, bool alloc = true)
        : KokkosVector<Particle<Model>, MemorySpace>(
              capacity, alloc, "particle_list")
    {
    }

    auto static with_capacity(std::size_t capacity)
    {
      auto rhs = ParticleList(capacity, false);
      return rhs;
    }

    KOKKOS_INLINE_FUNCTION Particle<Model> &operator[](size_t i)
    {
      return this->_owned_data(i);
    }

    KOKKOS_INLINE_FUNCTION const Particle<Model> &operator[](size_t i) const
    {
      return this->_owned_data(i);
    }

    auto d() const
    {
      return this->ddata();
    }

    ParticleList() = default;

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

    void init(double weight)
    {
      auto local_data = this->_owned_data;
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

    void _spawn_alloc(size_t n, double weight)
    {
      ParticleList<MemorySpace, Model> newlist(n);
      newlist.init(weight);
      this->insert(newlist);
    }
    MC::KPRNG rng_instance;
  private:
    
  };
} // namespace MC

#endif