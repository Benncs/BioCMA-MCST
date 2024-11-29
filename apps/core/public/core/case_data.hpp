#ifndef __CASE_DATA__HPP__
#define __CASE_DATA__HPP__

#include <simulation/feed_descriptor.hpp>
#include <common/execinfo.hpp>
#include <core/simulation_parameters.hpp>
#include <memory>
#include <optional>
#include <simulation/simulation.hpp>
#include <simulation/transitionner.hpp>



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
    std::unique_ptr<Simulation::FlowMapTransitioner> transitioner;

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
  void exec(CaseData &&case_data);

  std::optional<Core::CaseData> load(const ExecInfo &exec, const UserControlParameters &&params,std::optional<Simulation::Feed::SimulationFeed> feed=std::nullopt);

} // namespace Core

#endif