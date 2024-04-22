#include "mc/unit.hpp"
#include "models/simple_model.hpp"
#include <any>
#include <post_process.hpp>

#include <iostream>

/*
  NUMBER OF 
*/



void post_process(Simulation::SimulationUnit &simulation)
{
  std::cout<<"----END---"<<std::endl;
  auto model = std::any_cast<std::shared_ptr<SimpleModel>
  &>(simulation.mc_container->to_process[0].data);
  auto total_events = MC::EventContainer::reduce(simulation.mc_unit->ts_events);
  std::cout<<"-------"<<std::endl;
  std::cout<<"Death events: "<<total_events.get<MC::EventType::Death>()<<std::endl;
  std::cout<<"Division events: "<<total_events.get<MC::EventType::NewParticle>()<<std::endl;
  std::cout<<"mass: "<<model->xi.mass<<std::endl;
  std::cout<<"-------"<<std::endl;
  std::cout << simulation.mc_container->to_process.size() << '\n';

  auto d = simulation.mc_unit->domain.getDistribution();

  size_t count = 0;
  for (auto &&i : d)
  {
    std::cout << i << " ";
    count += i;
  }
  std::cout << '\n';
  std::cout << count << '\n';
  // TODO
}