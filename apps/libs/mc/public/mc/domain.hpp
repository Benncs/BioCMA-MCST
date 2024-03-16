#ifndef __MC_REACTORDOMAIN_HPP__
#define __MC_REACTORDOMAIN_HPP__

#include <common/common_types.hpp>
#include <mc/container_state.hpp>
#include <vector>

namespace MC
{
  class ReactorDomain
  {
  public:
    ReactorDomain() = default;
    ReactorDomain(ReactorDomain &&other);
    ReactorDomain(NumberView volumes,
                  std::vector<std::vector<size_t>> &&_neighbors);
    ReactorDomain(const ReactorDomain &other) = delete;

    inline auto &operator[](size_t i_c)
    {
      return this->containers[i_c];
    }
    inline auto &operator[](size_t i_c) const
    {
      return this->containers[i_c];
    }

    ReactorDomain &operator=(const ReactorDomain &other) = delete;

    double total_volume() const
    {
      return this->_total_volume;
    }
    ReactorDomain &operator=(ReactorDomain &&other) noexcept;

    size_t id;

  private:
    double _total_volume;

    std::vector<ContainerState> containers;
    std::vector<std::vector<size_t>> neighbors;
  };
} // namespace MC

#endif //__MC_REACTORDOMAIN_HPP__