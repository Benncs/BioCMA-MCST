#include "cma_read/neighbors.hpp"
#include "mc/particles/particles_list.hpp"
#include "mc/prng/prng.hpp"
#include <Kokkos_Core.hpp>
#include <Kokkos_Macros.hpp>
#include <Kokkos_Random.hpp>
#include <cstddef>
#include <decl/Kokkos_Declare_OPENMP.hpp>
#include <mc/mcinit.hpp>
#include <omp.h>
#include <utility>

namespace MC
{

  static double get_initial_weight(double scale_factor,
                                   double particle_concentration,
                                   double total_volume,
                                   size_t n_particles);

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
       DistributionVariantInt &&param,double x0)
  {

    auto unit = std::make_unique<MonteCarloUnit>();
    unit->domain = ReactorDomain(volumes, neighbors);

    size_t particle_per_process = n_particles / info.n_rank;

    size_t remainder = n_particles % info.n_rank;
    if (remainder != 0 && info.current_rank == info.n_rank - 1)
    {
      particle_per_process += remainder;
    }


    double weight =
        get_initial_weight(1., x0, unit->domain.getTotalVolume(), n_particles);

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

    // auto distribution =
    // MC::get_distribution_int<size_t>(distribution_variant);

    auto view_particle = p_list.data_span();
    Kokkos::View<MC::ContainerState *, Kokkos::LayoutRight> domain_view(
        domain.data().data(), n_compartment);

    const Kokkos::RangePolicy<> range(0, p_list.size());

    size_t max_compartment = (int)n_compartment;
    auto rng = std::make_shared<MC::KPRNG>();

    Kokkos::parallel_for(
        range, KOKKOS_LAMBDA(auto &&i_p) {
          auto &&particle = view_particle[i_p];

          particle.current_container = rng->uniform_u(
              0, max_compartment); // CAREFULL Kokkos rng in [0,n[
          particle.rng = rng;

          auto &i_container = domain_view[particle.current_container];
          init_kernel(particle);
          __ATOM_INCR__(i_container.n_cells);
        });
    Kokkos::fence();
  }

  double get_initial_weight(double scale_factor,
                            double particle_concentration,
                            double total_volume,
                            size_t n_particles)
  {
    //Scale factor is a fine tunning adjustement in case of misprediction of particle weight
    //particle_concentration is expected to be the real cell concentration in g/L (kg/m3)
    //Total volume is expected to be in m3
    //As a result we can calculate the mass carried by each MC particle 
    return scale_factor*particle_concentration*total_volume / static_cast<double>(n_particles);
  }

} // namespace MC