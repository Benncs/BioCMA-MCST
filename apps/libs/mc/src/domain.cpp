#include "mc/container_state.hpp"
#include <algorithm>
#include <cstddef>
#include <iterator>
#include <mc/domain.hpp>
#include <ranges>
#include <stdexcept>

//TODO REMOVE
#include <iostream>

namespace MC
{
  ReactorDomain::ReactorDomain(ReactorDomain &&other) noexcept
  {
    if (this != &other)
    {
      this->containers = std::move(other.containers);
      this->id = other.id;
      this->_total_volume = other._total_volume;
    }
  }

  void ReactorDomain::setVolumes(std::span<double const> volumesgas,
                                 std::span<double const> volumesliq)
  {
    // #pragma omp parallel for
    this->_total_volume = 0;
    for (size_t i_c = 0; i_c < volumesgas.size(); ++i_c)
    {
      containers[i_c].volume_liq = volumesliq[i_c];
      containers[i_c].volume_gas = volumesgas[i_c];
      this->_total_volume += volumesliq[i_c];
    }
  }
  ReactorDomain::ReactorDomain(std::span<double> volumes,
                               std::vector<std::vector<size_t>> &&_neighbors)
      :  neighbors(std::move(_neighbors))
  {
    double totv = 0.;
    std::transform(volumes.begin(),
                   volumes.end(),
                   std::back_inserter(this->containers),
                   [&totv, i = 0](auto &&v) mutable
                   {
                     auto c = ContainerState();
                     c.volume_liq = v;
                     c.id = i++;
                     totv += v;
                     return c;
                   });

    this->_total_volume = totv;
  }

  ReactorDomain &ReactorDomain::operator=(ReactorDomain &&other) noexcept
  {
    if (this != &other)
    {
      this->id = other.id;
      this->containers = std::move(other.containers);
      this->neighbors = std::move(other.neighbors);
      this->_total_volume = other._total_volume;
    }
    return *this;
  }

  std::vector<size_t> ReactorDomain::getDistribution() const
  {

    auto view =
        containers | std::views::transform([](const MC::ContainerState &cs)
                                           { return cs.n_cells; });

    return {view.begin(), view.end()};
  }

  ReactorDomain ReactorDomain::reduce(std::span<size_t> data,
                                      size_t original_size,
                                      size_t n_rank)
  {
    ReactorDomain reduced;
    reduced.containers.resize(original_size);

    if (data.size() != original_size * n_rank)
    {
      throw std::runtime_error("Cannot reduce different reactor type");
    }
    for (size_t i_rank = 0; i_rank < n_rank; ++i_rank)
    {
      for (size_t i_c = 0; i_c < original_size; ++i_c)
      {
        reduced.containers[i_c].n_cells += data[i_c + i_rank * original_size];
      }
    }

    return reduced;
  }

  // ReactorDomain ReactorDomain::reduce(std::span<ReactorDomain> data)
  // {
  //   ReactorDomain reduced;
  //   size_t size = data[0].containers.size();
  //   auto& containers = reduced.containers;
  //   for(auto&& domain : data)
  //   {
  //     if(size!=domain.n_compartments())
  //     {
  //       throw std::runtime_error("Cannot reduce diffferent reactor type");
  //     }
  //     for(size_t i =0 ;i<size;++i)
  //     {
  //       containers[i].n_cells+=domain.containers[i].n_cells;
  //     }
  //   }

  //   return reduced;
  // }

} // namespace MC