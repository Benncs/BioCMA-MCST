#include <cma_utils/discontinuous_fmt.hpp>
#include <cma_utils/transitionner.hpp>
namespace CmaUtils
{

  CmaUtils::PreCalculatedHydroState& DiscontinuousFMT::current_liq_hydro_state() noexcept
  {
    return liquid_pc[getFlowIndex()];
  }

  CmaUtils::PreCalculatedHydroState& DiscontinuousFMT::current_gas_hydro_state() noexcept
  {
    return gas_pc[getFlowIndex()];
  }

  [[nodiscard]] const CmaRead::ReactorState&
  DiscontinuousFMT::get_current_reactor_state() const noexcept
  {
    return get_unchecked(getFlowIndex());
  }

  DiscontinuousFMT::DiscontinuousFMT(std::size_t _n_flowmap,
                                     std::size_t _n_per_flowmap,
                                     std::size_t number_time_step,
                                     std::unique_ptr<CmaRead::FlowIterator>&& _iterator,
                                     bool is_two_phase_flow)
      : FlowMapTransitionner(
            _n_flowmap, _n_per_flowmap, number_time_step, std::move(_iterator), is_two_phase_flow)
  {
  }

  void DiscontinuousFMT::update_flow()
  {
    if (need_liquid_state())
    {
      calculate_full_state(
          get_current_reactor_state(), current_liq_hydro_state(), current_gas_hydro_state());
    }
  }

} // namespace CmaUtils