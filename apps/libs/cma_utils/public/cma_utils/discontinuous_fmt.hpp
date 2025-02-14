#ifndef __CMA_UTILS_D_FMT_HPP__
#define __CMA_UTILS_D_FMT_HPP__

#include <cma_utils/transitionner.hpp>

namespace CmaUtils
{
  class DiscontinuousFMT : public FlowMapTransitionner
  {
  public:
    DiscontinuousFMT(std::size_t _n_flowmap,
                     std::size_t _n_per_flowmap,
                     std::size_t number_time_step,
                     std::unique_ptr<CmaRead::FlowIterator>&& _iterator,
                     bool is_two_phase_flow);
    void update_flow() override;

  protected:
    [[nodiscard]] const CmaRead::ReactorState& get_current_reactor_state() const noexcept override;
    CmaUtils::PreCalculatedHydroState& current_liq_hydro_state() noexcept override;
    CmaUtils::PreCalculatedHydroState& current_gas_hydro_state() noexcept override;
  };
} // namespace CmaUtils

#endif