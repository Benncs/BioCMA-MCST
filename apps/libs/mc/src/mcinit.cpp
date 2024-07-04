#include "cma_read/neighbors.hpp"
#include "mc/particles/particles_list.hpp"
#include <mc/mcinit.hpp>
#include <omp.h>
#include <utility>

namespace MC
{

  void post_init_container(const std::function<void(MC::Particles &)>& init_kernel,
                           size_t n_thread,
                           ReactorDomain &domain,
                           ParticlesList &p_list,
                           size_t n_compartment,
                           DistributionVariantInt distribution_variant);

  std::unique_ptr<MonteCarloUnit>
  init(std::function<void(MC::Particles &)> init_kernel,
       const ExecInfo &info,
       size_t n_particles,
       std::span<double> volumes,
       CmaRead::Neighbors::Neighbors_const_view_t &neighbors,
       DistributionVariantInt &&param)
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
    // double weight = 1 / static_cast<double>(n_particles); // DUMMY INIT
    double x0 = 0.05; // g/l
    double v = 0.02;  // m3
    double m_part = 1e-15;

    double weight = (x0 * v) / (n_particles * m_part);
    // std::cout<<"we"<<weight<<std::endl;
    unit->container = ParticlesContainer(particle_per_process, weight);
    unit->extras.resize(info.thread_per_process);
    post_init_container(init_kernel,
                        info.thread_per_process,
                        unit->domain,
                        unit->container.to_process,
                        unit->domain.getNumberCompartments(),
                        std::move(param));
    return unit;
  }

  void post_init_container(const std::function<void(MC::Particles &)>& init_kernel,
                           size_t n_thread,
                           ReactorDomain &domain,
                           ParticlesList &p_list,
                           size_t n_compartment,
                           DistributionVariantInt distribution_variant)
  {

#pragma omp parallel default(none),                                            \
    shared(p_list, n_compartment, distribution_variant, domain,init_kernel),               \
    num_threads(n_thread)
    {
      auto prng = MC::PRNG::get_rng(omp_get_thread_num());
      auto distribution =
          MC::get_distribution_int<size_t>(distribution_variant);
      const auto size_p = p_list.size();
#pragma omp for schedule(static)
      for (size_t i_p = 0; i_p < size_p; i_p++)
      {

        auto &&particle = p_list[i_p]; //*it;
        particle.current_container = distribution(prng);

        auto &i_container = domain[particle.current_container];

        init_kernel(particle);
        __ATOM_INCR__(i_container.n_cells);
      }
    }
  }

} // namespace MC