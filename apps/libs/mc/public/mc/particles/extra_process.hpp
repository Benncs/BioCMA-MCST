#ifndef __MC__EXTRA_PROCESS_CONTAINER_HPP__
#define __MC__EXTRA_PROCESS_CONTAINER_HPP__

#include <mc/particles/mcparticles.hpp>
#include <mc/particles/particle_list.hpp>

namespace MC
{

  /**
   * @brief Helper class to handle data updated during kernel MonteCarlo process
   * @tparam MemorySpace: In which memory space the results are updated
   */
  template <typename MemorySpace, ParticleModel T> struct Results
  {

    using ResultsVIew =
        Kokkos::View<Results<MemorySpace, T>,
                     MemorySpace>; ///< Alias only for readibility

    MC::ParticleList<MemorySpace, T>
        extra_process; ///< Fixed size buffer with particles spawnt during the
                       ///< current cycle that will be process at the end of the
                       ///< kernel
    int waiting_allocation_particle =
        0; ///< Counter to know how many particle to allocate at the end at the
           ///< kernel if extra_process buffer is full.

    /**
     * @brief Main constructor
     * @param capacity: first size of buffer, must be adapated to expected
     * particle division number during cycle.
     * @warning According to implementation, buffer itself can allocate more
     * than capacity
     */
    explicit Results(size_t capacity)
        : extra_process(
              MC::ParticleList<MemorySpace, T>::with_capacity(capacity))
    {
    }

    Results() = default;

    /**
     * @brief Clear data after cycle to prepare next
     * @param new_capacity for new size of buffer, buffer doesn't shrink,
     * only grows
     */
    void clear(size_t new_capacity)
    {
      // TODO Handle case where user want to reduce capacity
      if (extra_process.size() < extra_process.capacity())
      {
        extra_process.resize(new_capacity);
      }
      extra_process.clear();
      waiting_allocation_particle = 0;
    }

    /**
     * @brief Update kernel view performing deepcopy
     */
    void update_view(ResultsVIew results_view)
    {
      Kokkos::deep_copy(results_view, *this);
    }

    /**
     * @brief Get correct view to be pass to kernel 
     */
    auto get_view()
    {
      Kokkos::View<Results<MemorySpace, T>, MemorySpace> results(
          "Results_view");
      update_view(results);
      return results;
    }


    // void from_view(Kokkos::View<Results<MemorySpace, T>, MemorySpace> _r_v)
    // {
    //   extra_process = _r_v().extra_proces;
    //   waiting_allocation_particle = _r_v().waiting_allocation_particle;
    // }
  };

} // namespace MC
#endif