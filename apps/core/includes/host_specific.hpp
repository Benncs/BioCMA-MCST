#ifndef __HOST_SPECIFIC_HPP__
#define __HOST_SPECIFIC_HPP__

// #include <common/common.hpp>

#include <common/logger.hpp>
#include <core/simulation_parameters.hpp>
#include <dataexporter/partial_exporter.hpp>
#include <memory>

// Foward declaration
namespace Simulation
{
  class SimulationUnit;
} // namespace Simulation

namespace CmaUtils
{
  class FlowMapTransitionner;
}

/**
 * @brief Main program executed on rank 0.
 *
 * This function is the main processing routine that is executed exclusively on
 * rank 0.
 * @param exec The `ExecInfo` object containing details about the execution
 * environment.
 * @param simulation The `Simulation::SimulationUnit` object representing the
 * simulation.
 * @param params The `SimulationParameters` object containing settings for the
 * simulation.
 * @param transitioner A unique pointer to the `Simulation::FlowMapTransitioner`
 * for handling flow map transitions.
 */
void host_process(std::shared_ptr<IO::Logger> logger, /*logger is passed by value because it can be
                                                         set into maindataexporter*/
                  const ExecInfo& exec,
                  Simulation::SimulationUnit& simulation,
                  const Core::SimulationParameters& params,
                  std::unique_ptr<CmaUtils::FlowMapTransitionner>&& transitioner,
                  Core::PartialExporter& partial_exporter);

#endif //__HOST_SPECIFIC_HPP__