#include <mc/mcinit.hpp>

namespace MC
{

  std::unique_ptr<ParticlesContainer> init_container(ExecInfo &info,
                                                     size_t n_particles)
  {
    size_t particle_per_process = n_particles / info.n_rank;

    size_t remainder = n_particles % info.n_rank;
    if (remainder != 0 && info.current_rank == info.n_rank - 1)
    {
      particle_per_process += remainder;
    }

    return std::make_unique<ParticlesContainer>(particle_per_process);
  }

  std::unique_ptr<MonteCarloUnit>
  init_unit(ExecInfo &info,
            NumberView volumes,
            std::vector<std::vector<size_t>> &&neighbors)
  {

    auto unit = std::make_unique<MonteCarloUnit>();
    unit->domain = ReactorDomain(volumes, std::move(neighbors));

    return unit;
  }

} // namespace MC