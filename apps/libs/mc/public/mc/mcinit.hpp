#ifndef __MC_INIT_HPP__
#define __MC_INIT_HPP__

#include <common/common_types.hpp>
#include <common/execinfo.hpp>
#include <mc/unit.hpp>

#include <mc/particles/particles_container.hpp>
#include <memory>


namespace MC
{
  std::unique_ptr<ParticlesContainer> init_container(ExecInfo &info,
                                                     size_t n_particles);
  std::unique_ptr<MonteCarloUnit>
  init_unit(ExecInfo &info,
            NumberView volumes,
            std::vector<std::vector<size_t>> &&neighbors);
} // namespace MC

#endif //__MC_INIT_HPP__