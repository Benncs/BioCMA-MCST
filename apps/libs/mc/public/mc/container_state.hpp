#ifndef __MC_CONTAINER_STATE_HPP__
#define __MC_CONTAINER_STATE_HPP__
#include <cstddef>

namespace MC
{
  struct ContainerState
  {
    double volume;
    size_t n_cells;
    size_t id;
  };
} // namespace MC

#endif //__MC_CONTAINER_STATE_HPP__