#include <algorithm>
#include <cstddef>
#include <iterator>
#include <mc/domain.hpp>

namespace MC
{
  ReactorDomain::ReactorDomain(ReactorDomain &&other)
  {
    if (this != &other)
    {
      this->containers = std::move(other.containers);
    }
  }
  ReactorDomain::ReactorDomain(NumberView volumes,
                               std::vector<std::vector<size_t>> &&_neighbors)
  {
    double __total_volume = 0;
    std::transform(volumes.begin(),
                   volumes.end(),
                   std::back_inserter(this->containers),
                   [&__total_volume, i = 0](auto &&v) mutable
                   {
                     auto c = ContainerState();
                     c.volume = v;
                     c.id = i++;
                     __total_volume += v;
                     return c;
                   });

    this->_total_volume = __total_volume;

    this->neighbors = std::move(_neighbors);
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

    std::transform(containers.begin(),containers.end(),std::back_inserter(distribution),[](auto&& cs){
      return cs.n_cells;
    });
    return distribution;
  }

} // namespace MC