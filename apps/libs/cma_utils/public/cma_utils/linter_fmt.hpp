#ifndef __CMA_UTILS_LINTER_FMT_HPP__
#define __CMA_UTILS_LINTER_FMT_HPP__

#include <cma_utils/transitionner.hpp>
#include <stdexcept>

namespace CmaUtils
{
  class LinterFMT : public FlowMapTransitionner
  {
  public:
    LinterFMT(std::size_t _n_flowmap,
              std::size_t _n_per_flowmap,
              std::size_t number_time_step,
              std::unique_ptr<CmaRead::FlowIterator>&& _iterator,
              bool is_two_phase_flow)
    {
      throw std::runtime_error("LinterFMT not implemented yet");
    }
    void update_flow() override
    {
      throw std::runtime_error("LinterFMT not implemented yet");
    }

  protected:
    [[nodiscard]] const CmaRead::ReactorState& get_current_reactor_state() const noexcept override
    {
      throw std::runtime_error("LinterFMT not implemented yet");
    }
    CmaUtils::PreCalculatedHydroState& current_liq_hydro_state() override
    {
      throw std::runtime_error("LinterFMT not implemented yet");
    }
    CmaUtils::PreCalculatedHydroState& current_gas_hydro_state() override
    {
      throw std::runtime_error("LinterFMT not implemented yet");
    }
  };
} // namespace CmaUtils

#endif