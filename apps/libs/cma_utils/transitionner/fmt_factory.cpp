#include <common/logger.hpp>
#include <discontinuous_fmt.hpp>
#include <linter_fmt.hpp>
#include <memory>
#include <stdexcept>
#include <transitionner/transitioner_factory.hpp>

namespace CmaUtils
{
  std::unique_ptr<FlowMapTransitionner>
  get_transitioner(const std::shared_ptr<IO::Logger>& logger,
                   FlowmapTransitionMethod method,
                   std::size_t n_flowmap,
                   std::size_t n_per_flowmap,
                   std::size_t number_time_step,
                   std::unique_ptr<CmaRead::FlowIterator>&& iterator,
                   bool is_two_phase_flow)
  {
    switch (method)
    {
    case (CmaUtils::FlowmapTransitionMethod::Discontinuous):
    {
      if (logger)
      {
        logger->print("FlowMapTransitionner", "Using Discontinous transition");
      }
      return std::make_unique<DiscontinuousFMT>(n_flowmap,
                                                n_per_flowmap,
                                                number_time_step,
                                                std::move(iterator),
                                                is_two_phase_flow);
    }
    case (CmaUtils::FlowmapTransitionMethod::InterpolationFO):
    {
      if (logger)
      {
        logger->print("FlowMapTransitionner",
                      "Using linear interpolation transition");
      }
      return std::make_unique<LinterFMT>(n_flowmap,
                                         n_per_flowmap,
                                         number_time_step,
                                         std::move(iterator),
                                         is_two_phase_flow);
    }
    default:
    {
      throw std::invalid_argument("FlowmapTransitionMethod not correct");
      __builtin_unreachable(); // TODO use c++23 cross plateform unreachabl
    }
    }
  }
} // namespace CmaUtils