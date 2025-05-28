#include <linter_fmt.hpp>
#include <discontinuous_fmt.hpp>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <transitionner/transitioner_factory.hpp>

namespace CmaUtils
{
  std::unique_ptr<FlowMapTransitionner>
  get_transitioner(FlowmapTransitionMethod method,
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
      std::cout << "Using Discontinous transition\r\n";
      return std::make_unique<DiscontinuousFMT>(
          n_flowmap, n_per_flowmap, number_time_step, std::move(iterator), is_two_phase_flow);
    }
    case (CmaUtils::FlowmapTransitionMethod::InterpolationFO):
    {
      std::cout << "Using linear interpolation transition\r\n";
      return std::make_unique<LinterFMT>(
          n_flowmap, n_per_flowmap, number_time_step, std::move(iterator), is_two_phase_flow);
    }
    default:
    {
      throw std::invalid_argument("FlowmapTransitionMethod not correct");
      __builtin_unreachable(); // TODO use c++23 cross plateform unreachabl
    }
    }
  }
} // namespace CmaUtils