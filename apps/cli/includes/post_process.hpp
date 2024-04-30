#ifndef __HOST_POST_PROCESS_HPP__
#define __HOST_POST_PROCESS_HPP__

#include <simulation/simulation.hpp>


void save_results(ExecInfo& exec,SimulationParameters &params,
                  Simulation::SimulationUnit &simulation);


void post_process(ExecInfo& exec,SimulationParameters &params,
                  Simulation::SimulationUnit &simulation);


void show(Simulation::SimulationUnit &simulation);

#endif //__HOST_POST_PROCESS _HPP__