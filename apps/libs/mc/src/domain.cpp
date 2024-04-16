#include <algorithm>
#include <cstddef>
#include <iterator>
#include <mc/domain.hpp>

#include <iostream>

namespace MC
{
  ReactorDomain::ReactorDomain(ReactorDomain &&other) noexcept
  {
    if (this != &other)
    {
      this->containers = std::move(other.containers);
    }
  }

  void ReactorDomain::setVolumes(std::span<double> volumesgas,
                                 std::span<double> volumesliq)
  {
// #pragma omp parallel for
    for (size_t i_c = 0; i_c < volumesgas.size(); ++i_c)
    {
      containers[i_c].volume_liq = volumesliq[i_c];
      containers[i_c].volume_gas = volumesgas[i_c];
    }
  }
  ReactorDomain::ReactorDomain(NumberView volumes,
                               std::vector<std::vector<size_t>> &&_neighbors)
      : neighbors(std::move(_neighbors))
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

  std::vector<size_t> ReactorDomain::getDistribution()
  {
    std::vector<size_t> distribution;
    distribution.reserve(this->n_compartments());

    std::transform(containers.begin(),
                   containers.end(),
                   std::back_inserter(distribution),
                   [](auto &&cs) { return cs.n_cells; });
    return distribution;
  }

} // namespace MC