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

  [[deprecated]] std::unique_ptr<MonteCarloUnit>
  init(const ExecInfo &info,
       size_t n_particles,
       std::span<double> volumes,
       CmaRead::Neighbors::Neighbors_const_view_t &neighbors);

  std::unique_ptr<MonteCarloUnit>
  init(std::function<void(MC::Particles &)> init_kernel,
       const ExecInfo &info,
       size_t n_particles,
       std::span<double> volumes,
       CmaRead::Neighbors::Neighbors_const_view_t &neighbors,
       DistributionVariantInt &&param,
       double x0);

  // std::unique_ptr<ParticlesContainer> init_container(ExecInfo &info,
  //                                                    size_t n_particles);
  // std::unique_ptr<MonteCarloUnit>
  // init_unit(ExecInfo &info,
  //           std::span<double> volumes,
  //           Neighbors::Neighbors_const_view_t& neighbors);
} // namespace MC

#endif //__MC_INIT_HPP__