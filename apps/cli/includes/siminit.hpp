#ifndef __PROCESS_INIT_HPP__
#define __PROCESS_INIT_HPP__

#include <memory>

#include <cma_read/flow_iterator.hpp>
#include <common/common.hpp>
#include <simulation/simulation.hpp>

Simulation::SimulationUnit
init_simulation(ExecInfo &info,
                SimulationParameters &params,
                std::shared_ptr<FlowIterator> &_flow_handle,
                KModel model,MC::DistributionVariantInt&& initial_particle_distribution);

#endif //__PROCESS_INIT_HPP__