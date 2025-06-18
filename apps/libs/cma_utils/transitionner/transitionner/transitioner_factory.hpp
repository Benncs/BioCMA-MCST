#ifndef __CMA_UTILS_TRANSITIONNER_FACTORY_HPP__
#define __CMA_UTILS_TRANSITIONNER_FACTORY_HPP__

#include "common/logger.hpp"
#include <memory>
#include <transitionner/transitionner.hpp>

namespace CmaUtils
{

  /**
   * @brief Transition factory according to the given method.
   *
   * Creates and returns a unique pointer to a FlowMapTransitionner based on the provided method
   * and configuration parameters.
   * @param method The method to use for transition.
   * @param n_flowmap The number of flow maps.
   * @param n_per_flowmap The number of iteration per flow map.
   * @param number_time_step The number of time steps.
   * @param iterator Optional unique pointer to a FlowIterator.
   * @param is_two_phase_flow Boolean flag indicating whether the flow is two-phase.
   * @return Unique pointer to FlowMapTransitionner.
   */
  std::unique_ptr<FlowMapTransitionner>
  get_transitioner(const std::shared_ptr<IO::Logger>& logger,FlowmapTransitionMethod method,
                   std::size_t n_flowmap,
                   std::size_t n_per_flowmap,
                   std::size_t number_time_step,
                   std::unique_ptr<CmaRead::FlowIterator>&& iterator = nullptr,
                   bool is_two_phase_flow = false);
} // namespace CmaUtils

#endif