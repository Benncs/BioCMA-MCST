// #ifndef __PROCESS_INIT_HPP__
// #define __PROCESS_INIT_HPP__

// #include <common/common.hpp>
// #include <memory>
// #include <core/simulation_parameters.hpp>
// #include <wrap_init_model_selector.hpp>
// // Forward declaration
// namespace Simulation
// {
//   class SimulationUnit;
//   class FlowMapTransitioner;
// }; // namespace Simulation

// /**
//  * @brief Initializes a simulation unit with the given parameters and
//  * transitioner.
//  *
//  * This function creates and configures a `Simulation::SimulationUnit` instance
//  * using the provided execution information, simulation parameters, and flow map
//  * transitioner.
//  *
//  * @param info The `ExecInfo` object containing details about the execution
//  * environment.
//  * @param params The `SimulationParameters` object with settings for the
//  * simulation.
//  * @param transitioner A unique pointer used to manage flowmap transitions.
//  * @return A unique pointer to the initialized `Simulation::SimulationUnit`.
//  */
// std::unique_ptr<Simulation::SimulationUnit>
// init_simulation(const ExecInfo &info,
//                 Core::SimulationParameters &params,
//                 std::unique_ptr<Simulation::FlowMapTransitioner> &transitioner);

// #endif //__PROCESS_INIT_HPP__