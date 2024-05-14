#include "cma_read/neighbors.hpp"
#include <mc/mcinit.hpp>
namespace MC
{

  std::unique_ptr<MonteCarloUnit> init(ExecInfo &info,
                      size_t n_particles,
                      std::span<double> volumes,
                      Neighbors::Neighbors_const_view_t &neighbors)
  { 

    auto unit = std::make_unique<MonteCarloUnit>();
    unit->domain = ReactorDomain(volumes, neighbors);
    unit->ts_events.resize(info.thread_per_process);

    size_t particle_per_process = n_particles / info.n_rank;

    size_t remainder = n_particles % info.n_rank;
    if (remainder != 0 && info.current_rank == info.n_rank - 1)
    {
      particle_per_process += remainder;
    }
    double weight = 1 / static_cast<double>(n_particles); // DUMMY INIT


    unit->container = ParticlesContainer(particle_per_process, 1 / weight);
    unit->extras.resize(info.thread_per_process);
    return unit;

  }

  // std::unique_ptr<ParticlesContainer> init_container(ExecInfo &info,
  //                                                    size_t n_particles)
  // {
  //   size_t particle_per_process = n_particles / info.n_rank;

  //   size_t remainder = n_particles % info.n_rank;
  //   if (remainder != 0 && info.current_rank == info.n_rank - 1)
  //   {
  //     particle_per_process += remainder;
  //   }

  //   // size_t total_population = 10'000'000; // FIXME
  //   // double weight = static_cast<double>(n_particles) /
  //   //                 static_cast<double>(total_population);

  //   // double weight = 1/(0.5*(90./1000.)/(double(n_particles)*1e-13));

  //   double weight = 1 / static_cast<double>(n_particles); // DUMMY INIT
  //   return std::make_unique<ParticlesContainer>(
  //       particle_per_process, 1 / weight);
  // }

  // std::unique_ptr<MonteCarloUnit>
  // init_unit(ExecInfo &info,
  //           std::span<double> volumes,
  //           Neighbors::Neighbors_const_view_t &neighbors)
  // {

  //   auto unit = std::make_unique<MonteCarloUnit>();
  //   unit->domain = ReactorDomain(volumes, neighbors);
  //   unit->ts_events.resize(info.thread_per_process);

  //   return unit;
  // }

} // namespace MC