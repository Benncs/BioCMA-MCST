#ifndef __MC_CONTAINER_STATE_HPP__
#define __MC_CONTAINER_STATE_HPP__
#include <cstddef>
#include <span> 


namespace MC
{
  struct ContainerState
  {
    double volume;
    size_t n_cells;
    size_t id;

    std::span<double const> concentrations;
  };
} // namespace MC

#endif //__MC_CONTAINER_STATE_HPP__