#ifndef __SYNC_NODES_HPP__
#define __SYNC_NODES_HPP__

#include <common/execinfo.hpp>
#include <simulation/simulation.hpp>

/**
 * @brief Synchronization after particle processing.
 *
 * This function performs synchronization operations for the simulation. It
 * ensures every computing unit is synchronized before going forward in time
 * step
 *
 * @param exec The `ExecInfo` object containing details about the execution
 * environment.
 * @param simulation The `Simulation::SimulationUnit` object representing the
 * simulation to be synchronized.
 */
void sync_step(const ExecInfo &exec, Simulation::SimulationUnit &simulation);

/**
 * @brief Synchronizes and resets simulation state for the next time step.
 *
 * This function performs necessary synchronization operations and resets
 * step-specific data, such as contributions and counters, in preparation for
 * the next time step in the simulation.
 *
 * @param exec The `ExecInfo` object containing details about the execution
 * environment.
 * @param simulation The `Simulation::SimulationUnit` object representing the
 * simulation being synchronized.
 */
void sync_prepare_next(Simulation::SimulationUnit &simulation);

/**
 * @brief Final synchronization before exporting results.
 *
 * This function performs the last synchronization and gathering operations
 * across computing units to merge results before exporting them.
 *
 * @param exec The `ExecInfo` object containing details about the execution
 * environment.
 * @param simulation The `Simulation::SimulationUnit` object representing the
 * simulation being synchronized.
 */
void last_sync(const ExecInfo &exec, Simulation::SimulationUnit &simulation);

#endif //__SYNC_NODES_HPP__