#include <discontinuous_fmt.hpp>
#include <linter_fmt.hpp>
#include <transitionner/proxy_cache.hpp>
#include <transitionner/transitionner.hpp>
namespace CmaUtils
{

  CmaUtils::ProxyPreCalculatedHydroState&
  DiscontinuousFMT::current_liq_hydro_state() noexcept
  {
    // With discontinuous we dont need to do anything just select the right
    // index
    return liquid_pc[getFlowIndex()];
  }

  CmaUtils::ProxyPreCalculatedHydroState&
  DiscontinuousFMT::current_gas_hydro_state() noexcept
  {
    // With discontinuous we dont need to do anything just select the right
    // index
    return gas_pc[getFlowIndex()];
  }

  [[nodiscard]] const CmaRead::ReactorState&
  DiscontinuousFMT::get_current_reactor_state() const noexcept
  {
    return get_unchecked(getFlowIndex());
  }

  DiscontinuousFMT::DiscontinuousFMT(
      std::size_t _n_flowmap,
      std::size_t _n_per_flowmap,
      std::size_t number_time_step,
      std::unique_ptr<CmaRead::FlowIterator>&& _iterator,
      bool is_two_phase_flow)
      : FlowMapTransitionner(_n_flowmap,
                             _n_per_flowmap,
                             number_time_step,
                             std::move(_iterator),
                             is_two_phase_flow)
  {
  }

  void DiscontinuousFMT::update_flow()
  {

    if (need_liquid_state_calculation())
    {
      // need_liquid_state_calculation ensures that this branch will be taken
      // only "n_flowmap" times
      calculate_full_state(get_current_reactor_state(),
                           current_liq_hydro_state(),
                           current_gas_hydro_state());
    }
  }

  CmaUtils::ProxyPreCalculatedHydroState&
  LinterFMT::current_liq_hydro_state() noexcept
  {
    throw std::runtime_error("LinterFMT not implemented yet");
  }
  CmaUtils::ProxyPreCalculatedHydroState&
  LinterFMT::current_gas_hydro_state() noexcept
  {
    throw std::runtime_error("LinterFMT not implemented yet");
  }

} // namespace CmaUtils