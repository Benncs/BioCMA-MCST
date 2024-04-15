#include "simulation/models/simple_model.hpp"
#include <any>
#include <post_process.hpp>

#include <iostream>

void post_process(Simulation::SimulationUnit &simulation)
{
  // auto model = std::any_cast<std::shared_ptr<SimpleModel> &>(simulation.container->to_process[0].data);
  // std::cout<<model->xi->mass<<std::endl;
  std::cout<<simulation.container->to_process.size()<<std::endl;

  auto d = simulation.mc_unit->domain.getDistribution();

  size_t count = 0;
  for(auto && i :d)
  {
    std::cout<<i<<" ";
    count+=i;
  }
  std::cout<<std::endl;
  std::cout<<count<<std::endl;
  //TODO
}