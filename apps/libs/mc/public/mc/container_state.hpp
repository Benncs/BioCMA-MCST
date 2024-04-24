#ifndef __MC_CONTAINER_STATE_HPP__
#define __MC_CONTAINER_STATE_HPP__
#include <cstddef>
#include <span>

namespace MC
{
  struct ContainerState
  {
    double volume_liq;
    double volume_gas;
    size_t n_cells; // TODO: BE CAREFUL WHEN DECREMENT
    size_t id;
    
    std::span<double const> concentrations;

  };
} // namespace MC

#endif //__MC_CONTAINER_STATE_HPP__