#ifndef __CMA_UTILS_TRANSITIONNER_FACTORY_HPP__
#define __CMA_UTILS_TRANSITIONNER_FACTORY_HPP__

#include <transitionner/transitionner.hpp>
#include <memory>

namespace CmaUtils
{
  
  

  std::unique_ptr<FlowMapTransitionner>
  get_transitioner(FlowmapTransitionMethod method,
                   std::size_t n_flowmap,
                   std::size_t n_per_flowmap,
                   std::size_t number_time_step,
                   std::unique_ptr<CmaRead::FlowIterator>&& iterator=nullptr,
                   bool is_two_phase_flow=false);
} // namespace CmaUtils

#endif