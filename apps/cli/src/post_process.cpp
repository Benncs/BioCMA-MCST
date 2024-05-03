#include "common/execinfo.hpp"
#include "mc/unit.hpp"
#include "models/light_model.hpp"
#include "models/simple_model.hpp"
#include <any>
#include <post_process.hpp>

#include <iostream>
#include <stdexcept>

/*
  NUMBER OF
*/

void model_specific(ExecInfo &exec,
                    SimulationParameters &params,
                    Simulation::SimulationUnit &simulation)
{
#ifndef USE_PYTHON_MODULE
  auto &model =
      std::any_cast<SimpleModel &>(simulation.mc_container->to_process[0].data);
  std::cout << "mass: " << model.xi.mass << std::endl;
#else
  // simulation.getModel().f_dbg(simulation.mc_container->to_process[0]);
#endif
}

void post_process(ExecInfo &exec,
                  SimulationParameters &params,
                  Simulation::SimulationUnit &simulation)
{

  const size_t process_size = simulation.mc_container->to_process.size();

  std::cout << "----END---" << std::endl;

  // FIXME
  try
  {
    model_specific(exec, params, simulation);
  }
  catch (...)
  {
    ///
  }

  auto total_events = simulation.mc_unit->ts_events[0];
  const size_t death_events = total_events.get<MC::EventType::Death>();
  const size_t new_events = total_events.get<MC::EventType::NewParticle>();

  auto d = simulation.mc_unit->domain.getDistribution();
  size_t count = 0;
  for (auto &&i : d)
  {
    std::cout << i << " ";
    count += i;
  }
  std::cout << '\n';

  std::cout << "\r\n-------\r\n";
  std::cout << "Death events: " << death_events << std::endl;
  std::cout << "Division events: " << new_events << std::endl;

  std::cout << "Starting number particle to process: " << params.n_particles
            << std::endl;
  std::cout << "Ending number particle to process: "
            << process_size * exec.n_rank << "(" << process_size << "*"
            << exec.n_rank << ")" << std::endl;

  std::cout << "Number living particle : " << count << std::endl;
  std::cout << "\r\n-------\r\n" << std::endl;
  ;
  if (exec.n_rank == 1)
  {
    assert(count == process_size - death_events);
    assert(count == params.n_particles + new_events - death_events);
  }

  // TODO
}

void show(Simulation::SimulationUnit &simulation)
{

  std::vector<double> mass(simulation.mc_unit->domain.getNumberCompartments());
  double totmass = 0.;

  auto d = simulation.mc_unit->domain.getDistribution();
  size_t count = 0;
  for (auto &&i : d)
  {
    std::cout << i << " ";
    count += i;
  }
  std::cout << '\n';

  // std::for_each(simulation.mc_container->to_process.begin(),
  //               simulation.mc_container->to_process.end(),
  //               [&mass, &totmass](auto &&p)
  //               {
  //                 auto &model = std::any_cast<LightModel &>(p.data);
  //                 totmass += model.mass * p.weight;
  //                 mass[p.current_container] += model.mass * p.weight;
  //               });

  std::cout << simulation.mc_unit->domain.getTotalVolume() * 1000 << std::endl;
  // auto concentration = totmass / (simulation.mc_unit->domain.getTotalVolume());
  // std::cout << "total mass: " << totmass << "\r\n"
  //           << "Mean concentration: " << concentration << "\r\n\r\\n";

  // for (size_t i = 0; i < simulation.mc_unit->domain.getNumberCompartments();
  //      ++i)
  // {
  //   double bio_concentrations =
  //       mass[i] / (simulation.mc_unit->domain[i].volume_liq);
  //   std::cout << bio_concentrations << " | ";
  // }
  // std::cout << std::endl;
}


void save_results(ExecInfo& exec,SimulationParameters &params,
                  Simulation::SimulationUnit &simulation)
                  {
                    
                  }