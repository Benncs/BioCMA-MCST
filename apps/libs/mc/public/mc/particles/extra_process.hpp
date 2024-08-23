#ifndef __MC__EXTRA_PROCESS_CONTAINER_HPP__
#define __MC__EXTRA_PROCESS_CONTAINER_HPP__

#include <mc/particles/mcparticles.hpp>
#include <mc/particles/particle_list.hpp>

namespace MC
{

  template <typename MemorySpace, MC::ParticleModel T> struct Results
  {

    using ResultsVIew = Kokkos::View<Results<MemorySpace, T>, MemorySpace>;

    MC::ParticleList<MemorySpace, T> extra_proces;
    int waiting_allocation_particle = 0;

    explicit Results(size_t capacity)
        : extra_proces(
              MC::ParticleList<MemorySpace, T>::with_capacity(capacity))
    {
    }

    Results() = default;

    void clear(size_t capacity)
    {
      extra_proces.resize(capacity);
      extra_proces.clear();
      waiting_allocation_particle = 0;
    }

    void update_view(ResultsVIew results_view)
    {
      Kokkos::deep_copy(results_view, *this);
    }

    auto get_view()
    {
      Kokkos::View<Results<MemorySpace, T>, MemorySpace> results(
          "Results_view");
      update_view(results);
      return results;
    }

    void from_view(Kokkos::View<Results<MemorySpace, T>, MemorySpace> _r_v)
    {
      extra_proces = _r_v().extra_proces;
      waiting_allocation_particle = _r_v().waiting_allocation_particle;
    }
  };

} // namespace MC
#endif