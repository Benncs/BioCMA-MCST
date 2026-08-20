#ifndef __SIMULATIONS_KERNEL_REDUCER_HPP__
#define __SIMULATIONS_KERNEL_REDUCER_HPP__

#include <Kokkos_Core.hpp>
#include <biocma_cst_config.hpp>
#include <cassert>
#include <common/common.hpp>
#include <mc/alias.hpp>
#include <mc/events.hpp>
#include <mc/traits.hpp>

namespace Simulation::KernelInline
{
  struct CycleReducerType
  {
    std::size_t waiting_allocation_particle;
    std::size_t dead_total;
    // std::size_t division;

    KOKKOS_INLINE_FUNCTION CycleReducerType&
    operator+=(const CycleReducerType& a)
    {
      this->waiting_allocation_particle += a.waiting_allocation_particle;
      this->dead_total += a.dead_total;
      //   this->division+=a.division;
      return *this;
    }
  };

  template <class Space> class CycleReducer
  {
  public:
    // Required for Concept
    using reducer = CycleReducer;
    using value_type = CycleReducerType;
    using result_view_type = Kokkos::View<value_type, Space>;

    KOKKOS_INLINE_FUNCTION
    void
    join(value_type& dest, const value_type& src) const
    {
      // dest.dead_total += src.dead_total;
      // dest.waiting_allocation_particle += src.waiting_allocation_particle;
      dest += src;
    }

    [[nodiscard]] KOKKOS_INLINE_FUNCTION value_type&
    reference() const
    {
      return *value.data();
    }

    KOKKOS_INLINE_FUNCTION
    result_view_type
    view() const
    {
      return value;
    }

    [[nodiscard]] KOKKOS_INLINE_FUNCTION bool
    references_scalar() const
    {
      return references_scalar_v;
    }

    // Optional
    KOKKOS_INLINE_FUNCTION
    void
    init(value_type& val) const
    {
      val.dead_total = 0;
      val.waiting_allocation_particle = 0;
      //   val.division=0;
    }

    // KOKKOS_INLINE_FUNCTION
    // void final(value_type& val) const
    // {
    //   // NOP
    // }

    // Part of Build-In reducers for Kokkos
    KOKKOS_INLINE_FUNCTION
    explicit CycleReducer(value_type& value_)
        : value(&value_), references_scalar_v(true)
    {
    }

    KOKKOS_INLINE_FUNCTION
    explicit CycleReducer(const result_view_type& value_)
        : value(value_), references_scalar_v(false)
    {
    }

  private:
    result_view_type value;
    bool references_scalar_v;
  };
}

#endif
