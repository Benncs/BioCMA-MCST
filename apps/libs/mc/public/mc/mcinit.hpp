#ifndef __MC_INIT_HPP__
#define __MC_INIT_HPP__

#include <common/common_types.hpp>
#include <common/execinfo.hpp>
#include <common/thread_safe_container.hpp>
#include <mc/particles/particles_container.hpp>
#include <mc/unit.hpp>
#include <memory>

namespace MC
{

  //   [[deprecated]] std::unique_ptr<MonteCarloUnit>
  //   init(const ExecInfo &info,
  //        size_t n_particles,
  //        std::span<double> volumes,
  //        CmaRead::Neighbors::Neighbors_const_view_t &neighbors);

  //   std::unique_ptr<MonteCarloUnit>
  //   init(std::function<void(MC::Particles &)> init_kernel,
  //        const ExecInfo &info,
  //        size_t n_particles,
  //        std::span<double> volumes,
  //        CmaRead::Neighbors::Neighbors_const_view_t &neighbors,
  //        DistributionVariantInt &&param,
  //        double x0);

  namespace
  {
    inline double get_initial_weight(double scale_factor,
                                     double particle_concentration,
                                     double total_volume,
                                     size_t n_particles)
    {
      // Scale factor is a fine tunning adjustement in case of misprediction of
      // particle weight particle_concentration is expected to be the real cell
      // concentration in g/L (kg/m3) Total volume is expected to be in m3 As a
      // result we can calculate the mass carried by each MC particle
      return scale_factor * particle_concentration * total_volume /
             static_cast<double>(n_particles);
    }
  } // namespace

  template <MC::ParticleModel Model>
  std::unique_ptr<MonteCarloUnit> inite(const ExecInfo &info,
                                        size_t n_particles,std::span<double> volumes,
                                        CmaRead::Neighbors::Neighbors_const_view_t neighbors,
                                        DistributionVariantInt &&param,
                                        double x0)
  {
    auto unit = std::make_unique<MonteCarloUnit>();

    unit->domain =
        ReactorDomain(volumes, neighbors);

    size_t particle_per_process = n_particles / info.n_rank;

    size_t remainder = n_particles % info.n_rank;
    if (remainder != 0 && info.current_rank == info.n_rank - 1)
    {
      particle_per_process += remainder;
    }

    double weight =
        get_initial_weight(1., x0, unit->domain.getTotalVolume(), n_particles);

    unit->container = ParticlesContainer<Model>(particle_per_process, weight);
    return unit;
  }

  // std::unique_ptr<ParticlesContainer> init_container(ExecInfo &info,
  //                                                    size_t n_particles);
  // std::unique_ptr<MonteCarloUnit>
  // init_unit(ExecInfo &info,
  //           std::span<double> volumes,
  //           Neighbors::Neighbors_const_view_t& neighbors);
} // namespace MC

#endif //__MC_INIT_HPP__