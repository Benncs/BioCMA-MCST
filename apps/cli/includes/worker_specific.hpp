#ifndef __WORKER_SPECIFIC_HPP__
#define __WORKER_SPECIFIC_HPP__


#include <common/execinfo.hpp>
#include <simulation/simulation.hpp>
#include <simulation/update_flows.hpp>

void workers_process(
    const ExecInfo &exec,
    Simulation::SimulationUnit &&simulation,
    const SimulationParameters &params,
    std::unique_ptr<Simulation::FlowMapTransitioner> &&transitioner);



#endif //__WORKER_SPECIFIC_HPP__