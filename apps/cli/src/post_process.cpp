#include "mc/unit.hpp"
#include "models/simple_model.hpp"
#include <any>
#include <post_process.hpp>

#include <iostream>
#include <stdexcept>

/*
  NUMBER OF
*/

void post_process(SimulationParameters &params,
                  Simulation::SimulationUnit &simulation)
{

  const size_t process_size = simulation.mc_container->to_process.size();

  std::cout << "----END---" << std::endl;

  std::cout << "-------" << std::endl;

  // FIXME
#ifndef USE_PYTHON_MODULE
  auto model = std::any_cast<std::shared_ptr<SimpleModel> &>(
      simulation.mc_container->to_process[0].data);
  std::cout << model->xi.mass << std::endl;
#else
  simulation.getModel().f_dbg(simulation.mc_container->to_process[0]);
#endif

  auto total_events =
      MC::EventContainer::reduce_local(simulation.mc_unit->ts_events);
  const size_t death_events = total_events.get<MC::EventType::Death>();
  const size_t new_events = total_events.get<MC::EventType::NewParticle>();

  std::cout << "Death events: " << death_events << std::endl;
  std::cout << "Division events: " << new_events << std::endl;

  std::cout << "-------" << std::endl;
  std::cout << process_size << '\n';

  auto d = simulation.mc_unit->domain.getDistribution();
  size_t count = 0;
  for (auto &&i : d)
  {
    std::cout << i << " ";
    count += i;
  }
  std::cout << '\n';
  std::cout << count << std::endl;
  assert(count == process_size - death_events);
  assert(count == params.n_particles + new_events - death_events);

  // TODO
}