#ifndef __HOST_POST_PROCESS_HPP__
#define __HOST_POST_PROCESS_HPP__

#include <simulation/simulation.hpp>
#include <data_exporter.hpp>

void save_results(ExecInfo& exec,SimulationParameters &params,
                  Simulation::SimulationUnit &simulation);


void post_process(ExecInfo& exec,SimulationParameters &params,
                  Simulation::SimulationUnit &simulation,std::unique_ptr<DataExporter>& exporter);


void show(Simulation::SimulationUnit &simulation);

#endif //__HOST_POST_PROCESS _HPP__