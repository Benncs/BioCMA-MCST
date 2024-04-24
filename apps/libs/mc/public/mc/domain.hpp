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
    size_t id{};

    ReactorDomain() = default;
    ReactorDomain(ReactorDomain &&other) noexcept;
    ReactorDomain(NumberView volumes,
                  std::vector<std::vector<size_t>> &&_neighbors);
    ReactorDomain(const ReactorDomain &other) = delete;
    ~ReactorDomain() = default;

    inline auto &operator[](size_t i_c)
    {
      return this->containers[i_c];
    }
    inline auto &operator[](size_t i_c) const
    {
      return this->containers[i_c];
    }

    ReactorDomain &operator=(const ReactorDomain &other) = delete;
    ReactorDomain &operator=(ReactorDomain &&other) noexcept;

    [[nodiscard]] decltype(auto) begin() const
    {
      return containers.begin();
    }
    [[nodiscard]] decltype(auto) end() const
    {
      return containers.end();
    }

    decltype(auto) begin()
    {
      return containers.begin();
    }
    decltype(auto) end()
    {
      return containers.end();
    }

    [[nodiscard]] inline auto n_compartments() const
    {
      return containers.size();
    }

    [[nodiscard]] inline double total_volume() const
    {
      return this->_total_volume;
    }

    void setVolumes(std::span<double const> volumesgas,
                    std::span<double const> volumesliq);

    void setLiquidNeighbors(const std::vector<std::vector<size_t>> &data)
    {
      neighbors = data;
    }

    [[nodiscard]] inline const std::vector<std::vector<size_t>> &
    getNeighbors() const
    {
      return neighbors;
    }

    [[nodiscard]] inline std::span<const size_t> getNeighbors(size_t i) const
    {
      return neighbors[i];
    }

    [[nodiscard]] std::vector<size_t> getDistribution() const;

  private:
    double _total_volume{};

    std::vector<ContainerState> containers;
    std::vector<std::vector<size_t>> neighbors;
  };

} // namespace MC

#endif //__MC_REACTORDOMAIN_HPP__