#ifndef __CMA_UTILS_D_FMT_HPP__
#define __CMA_UTILS_D_FMT_HPP__

#include <transitionner/transitionner.hpp>

namespace CmaUtils
{
  class DiscontinuousFMT final: public FlowMapTransitionner
  {
  public:
    DiscontinuousFMT(std::size_t _n_flowmap,
                     std::size_t _n_per_flowmap,
                     std::size_t number_time_step,
                     std::unique_ptr<CmaRead::FlowIterator>&& _iterator,
                     bool is_two_phase_flow);
    void update_flow() final;

  protected:
    [[nodiscard]] const CmaRead::ReactorState& get_current_reactor_state() const noexcept final;
    CmaUtils::ProxyPreCalculatedHydroState& current_liq_hydro_state() noexcept final;
    CmaUtils::ProxyPreCalculatedHydroState& current_gas_hydro_state() noexcept final;
  };
} // namespace CmaUtils

#endif