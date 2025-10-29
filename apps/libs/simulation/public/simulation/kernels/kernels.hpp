#ifndef __SIMULATION_KERNELS_HPP__
#define __SIMULATION_KERNELS_HPP__

#include <mc/unit.hpp>
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

    void launch_move(const std::size_t n_particle) const
    {

      if (move_kernel.enable_move)
      {
        const auto _policy_move =
            MC::get_policy<move_kernel_type, KernelInline::TagMove>(
                move_kernel, n_particle, false);

        Kokkos ::parallel_for("cycle_move", _policy_move, move_kernel);
      }
      if (move_kernel.enable_leave)
      {

        const auto _policy_leave =
            Kokkos::TeamPolicy<ComputeSpace, KernelInline::TagLeave>(
                32, Kokkos::AUTO);

        Kokkos ::parallel_reduce(
            "cycle_move_leave", _policy_leave, move_kernel, move_reducer);
      }
    }

    void launch_model(const std::size_t n_particle) const
    {
      constexpr bool is_reduce = true;
      const auto _policy =
          MC::get_policy<cycle_kernel_type, KernelInline::TagSecondPass>(
              cycle_kernel, n_particle, is_reduce);

      const auto scatter_policy =
          Kokkos::RangePolicy<KernelInline::TagFirstPass>(0, n_particle);

      Kokkos::parallel_reduce(
          "cycle_model",
          _policy,
          cycle_kernel,
          KernelInline::CycleReducer<ComputeSpace>(cycle_reducer));
      Kokkos::fence(); // TODO needed ?
      Kokkos::parallel_for("cycle_scatter", scatter_policy, cycle_kernel);
    }
  };

} // namespace Simulation::KernelInline

#endif
