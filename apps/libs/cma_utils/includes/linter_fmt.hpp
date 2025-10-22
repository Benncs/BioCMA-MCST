#ifndef __CMA_UTILS_LINTER_FMT_HPP__
#define __CMA_UTILS_LINTER_FMT_HPP__

#include <stdexcept>
#include <transitionner/transitionner.hpp>

namespace CmaUtils
{
  /**
 @brief Derived class to handle Linear interpolation transition between flowmaps

 Linear interpolation means that we pass don't pass from flowmap i to flormap
 i+1 but we use intermediate flowmap created with linear interpolation
  */
  class LinterFMT final : public FlowMapTransitionner
  {
  public:
    LinterFMT(std::size_t _n_flowmap,
              std::size_t _n_per_flowmap,
              std::size_t number_time_step,
              std::unique_ptr<CmaRead::FlowIterator>&& _iterator,
              bool is_two_phase_flow)
    {
      (void)_n_flowmap;
      (void)_n_per_flowmap;
      (void)number_time_step;
      (void)_iterator;
      (void)is_two_phase_flow;
      throw std::runtime_error("LinterFMT not implemented yet");
    }
    void update_flow() final
    {
      throw std::runtime_error("LinterFMT not implemented yet");
    }

  protected:
    [[nodiscard]] const CmaRead::ReactorState&
    get_current_reactor_state() const noexcept final
    {
      throw std::runtime_error("LinterFMT not implemented yet");
    }
    CmaUtils::ProxyPreCalculatedHydroState&
    current_liq_hydro_state() noexcept final;
    CmaUtils::ProxyPreCalculatedHydroState&
    current_gas_hydro_state() noexcept final;
  };
} // namespace CmaUtils

#endif