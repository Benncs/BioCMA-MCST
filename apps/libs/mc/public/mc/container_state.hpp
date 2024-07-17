#ifndef __MC_CONTAINER_STATE_HPP__
#define __MC_CONTAINER_STATE_HPP__
#include "common/execinfo.hpp"
#include <cstddef>
#include <cstdint>
#include <span>


namespace MC
{
  struct alignas(ExecInfo::cache_line_size) ContainerState
  {
    double volume_liq;
    double volume_gas;
    size_t n_cells; // TODO: BE CAREFUL WHEN DECREMENT
    uint64_t id;
    std::span<double const> concentrations;
  };
} // namespace MC

#endif //__MC_CONTAINER_STATE_HPP__