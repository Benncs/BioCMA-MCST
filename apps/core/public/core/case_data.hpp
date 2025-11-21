#ifndef __CASE_DATA__HPP__
#define __CASE_DATA__HPP__

#include "cma_utils/d_transitionner.hpp"
#include <common/execinfo.hpp>
#include <common/logger.hpp>
#include <core/simulation_parameters.hpp>
#include <memory>
#include <optional>
#include <simulation/feed_descriptor.hpp>
#include <simulation/mass_transfer.hpp>

namespace Simulation
{
  class SimulationUnit;
} // namespace Simulation

namespace CmaUtils
{
  class FlowMapTransitionner;
} // namespace CmaUtils

/**
  @brief Core component to perform simulation
 */
namespace Core
{

  /**
   * @brief Holds the data required to execute a simulation case.
   *
   * The `CaseData` struct encapsulates all the necessary components needed to
   * perform a simulation. It includes the simulation unit, parameters,
   * transitioner, and execution information. This structure is typically used
   * to manage and transfer simulation data between different stages or
   * components of the simulation process.
   */
  struct CaseData
  {

    CaseData();

    /**
     * @brief Unique pointer to the simulation unit.
     *
     * This is the core unit of the simulation, responsible for executing the
     * main simulation logic.
     */
    std::unique_ptr<Simulation::SimulationUnit> simulation;

    /**
     * @brief Parameters that configure the simulation.
     *
     * The simulation parameters control various aspects of the simulation,
     * such as time steps, boundary conditions, and solver configurations.
     */
    SimulationParameters params;

    /**
     * @brief Unique pointer to the flow map transitioner.
     *
     * The transitioner manages transitions between different flow states within
     * the simulation, ensuring consistency and accuracy in the simulation's
     * progression.
     */
    // std::unique_ptr<CmaUtils::FlowMapTransitionner> transitioner;

    CmaUtils::TransitionnerPtrType transitioner;

    /**
     * @brief Information about the execution environment.
     *
     * The execution information contains details such as the execution context,
     * environment settings, and other relevant metadata necessary for the
     * simulation run.
     */
    ExecInfo exec_info;
  };

  /**
   * @brief Start simulation
   */
  void exec(std::shared_ptr<IO::Logger> logger, CaseData&& case_data);

  std::optional<Core::CaseData>
  load(const ExecInfo& exec,
       const UserControlParameters&& params,
       std::optional<Simulation::Feed::SimulationFeed> feed = std::nullopt);

  /**
   * @brief Initializes the runtime environment based on command-line arguments
   * and simulation parameters.
   *
   * This function sets up the necessary runtime environment for the simulation
   * by:
   * - Initializing MPI (Message Passing Interface) if applicable.
   * - Setting up Kokkos for parallel programming.
   * - Configuring functions to be executed upon program exit.
   * - Handling signals to ensure proper shutdown and resource cleanup.
   *
   * The function uses the provided command-line arguments and simulation
   * parameters to configure the runtime environment accordingly.
   *
   * @param argc The number of command-line arguments.
   * @param argv The array of command-line arguments.
   * settings for the simulation.
   * @return An `ExecInfo` object containing details about the initialized
   * runtime environment, including execution context and other relevant
   * metadata.
   */
  ExecInfo runtime_init(int argc,
                        char** argv,
                        std::optional<std::size_t> force_run_id = std::nullopt);
} // namespace Core

#endif //!__CASE_DATA__HPP__
