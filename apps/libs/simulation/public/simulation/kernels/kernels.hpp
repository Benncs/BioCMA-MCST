#ifndef __SIMULATION_KERNELS_HPP__
#define __SIMULATION_KERNELS_HPP__

#include <simulation/kernels/model_kernel.hpp>
#include <simulation/kernels/move_kernel.hpp>

namespace Simulation::KernelInline
{

  template <typename Space, ModelType Model> struct Functors
  {
    using FModel = Model;

    using cycle_reducer_view_type =
        KernelInline::CycleReducer<Space>::result_view_type;
    using move_reducer_view_type = Kokkos::View<std::size_t, Space>;
    using cycle_kernel_type = CycleFunctor<Model>;
    using move_kernel_type = MoveFunctor;

    cycle_reducer_view_type cycle_reducer;
    move_reducer_view_type move_reducer;
    cycle_kernel_type cycle_kernel;
    move_kernel_type move_kernel;

    explicit Functors(cycle_kernel_type&& ck, move_kernel_type&& mk)
        : cycle_reducer("cycle_reducer"), move_reducer("move_reducer"),
          cycle_kernel(std::move(ck)), move_kernel(std::move(mk))
    {
    }
  };
} // namespace Simulation::KernelInline

#endif
