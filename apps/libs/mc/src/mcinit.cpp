#include "cma_read/neighbors.hpp"
#include "mc/particles/particles_list.hpp"
#include <Kokkos_Core.hpp>
#include <Kokkos_Macros.hpp>
#include <Kokkos_Random.hpp>
#include <decl/Kokkos_Declare_OPENMP.hpp>
#include <mc/mcinit.hpp>
#include <omp.h>
#include <utility>

namespace MC
{

  void
  post_init_container(const std::function<void(MC::Particles &)> &init_kernel,
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

    size_t particle_per_process = n_particles / info.n_rank;

    size_t remainder = n_particles % info.n_rank;
    if (remainder != 0 && info.current_rank == info.n_rank - 1)
    {
      particle_per_process += remainder;
    }

    double x0 = 5;                            // g/l
    double v = unit->domain.getTotalVolume(); // m3
    double m_part = 1e-15;
    // double weight = (x0 * v) / (n_particles); // (x0 * v) / (n_particles *
    // m_part); // double weight = (x0 *
    // v)/static_cast<double>(n_particles);//(x0 * v)/n_particles;
    double weight = 1 / static_cast<double>(n_particles);

    unit->container = ParticlesContainer(particle_per_process, weight);
    // unit->extras.resize(info.thread_per_process);
    post_init_container(init_kernel,
                        info.thread_per_process,
                        unit->domain,
                        unit->container.to_process,
                        unit->domain.getNumberCompartments(),
                        param);
    return unit;
  }

  void
  post_init_container(const std::function<void(MC::Particles &)> &init_kernel,
                      size_t n_thread,
                      ReactorDomain &domain,
                      ParticlesList &p_list,
                      size_t n_compartment,
                      DistributionVariantInt distribution_variant)
  {

    auto distribution = MC::get_distribution_int<size_t>(distribution_variant);

    auto view_particle = p_list.data_span();
    Kokkos::View<MC::ContainerState *, Kokkos::LayoutRight> domain_view(
        domain.data().data(), domain.getNumberCompartments());

    const Kokkos::RangePolicy<> range(0, p_list.size());

    auto max_compartment = domain.getNumberCompartments() - 1;

    Kokkos::Random_XorShift1024_Pool<> random_pool(12345);

    Kokkos::parallel_for(
        range, KOKKOS_LAMBDA(auto &&i_p) {
          auto seed = std::random_device{}();
          auto gen = std::mt19937(seed);
          auto &&particle = view_particle[i_p];
          // particle.current_container = distribution(gen);
          auto generator = random_pool.get_state();

          particle.current_container = generator.urand(0, max_compartment);

          random_pool.free_state(generator);
          auto &i_container = domain_view[particle.current_container];
          init_kernel(particle);
          __ATOM_INCR__(i_container.n_cells);
        });
    Kokkos::fence();

    // std::vector<MC::PRNG> rng(n_thread);
    // #pragma omp parallel default(none), \
    // shared(p_list, \
    //            n_compartment, \
    //            distribution_variant, \
    //            domain, \
    //            init_kernel, \
    //            rng),                                                       \
    // num_threads(n_thread)
    //     {
    //       auto &prng = rng[omp_get_thread_num()].rng();
    //       auto distribution =
    //           MC::get_distribution_int<size_t>(distribution_variant);
    //       const auto size_p = p_list.size();
    // #pragma omp for schedule(static)
    //       for (size_t i_p = 0; i_p < size_p; i_p++)
    //       {

    //         auto &&particle = p_list[i_p]; //*it;
    //         particle.current_container = distribution(prng);

    //         auto &i_container = domain[particle.current_container];

    //         init_kernel(particle);
    //         // __ATOM_INCR__(i_container.n_cells);
    //         #pragma omp atomic
    //         i_container.n_cells++;
    //       }
    //     }
  }

} // namespace MC