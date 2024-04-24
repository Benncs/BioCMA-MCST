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

    [[nodiscard]] decltype(auto) begin() const;
    [[nodiscard]] decltype(auto) end() const;

    decltype(auto) begin();
    decltype(auto) end();

    [[nodiscard]] auto n_compartments() const;

    [[nodiscard]] double total_volume() const;

    void setVolumes(std::span<double const> volumesgas,
                    std::span<double const> volumesliq);

    void setLiquidNeighbors(const std::vector<std::vector<size_t>> &data);

    [[nodiscard]] const std::vector<std::vector<size_t>> &getNeighbors() const;

    [[nodiscard]] std::span<const size_t> getNeighbors(size_t i) const;

    [[nodiscard]] std::vector<size_t> getDistribution() const;

    static ReactorDomain reduce(std::span<size_t> data,size_t original_size,size_t n_rank);

    std::span<ContainerState> data();

  private:
    double _total_volume{};

    std::vector<ContainerState> containers;
    std::vector<std::vector<size_t>> neighbors;
  };

  inline std::span<ContainerState> ReactorDomain::data()
  {
    return containers;
  }

  inline std::span<const size_t> ReactorDomain::getNeighbors(size_t i) const
  {
    return neighbors[i];
  }

  inline const std::vector<std::vector<size_t>> &
  ReactorDomain::getNeighbors() const
  {
    return neighbors;
  }

  inline void ReactorDomain::setLiquidNeighbors(
      const std::vector<std::vector<size_t>> &data)
  {
    neighbors = data;
  }

  inline auto ReactorDomain::n_compartments() const
  {
    return containers.size();
  }

  inline double ReactorDomain::total_volume() const
  {
    return this->_total_volume;
  }

  inline decltype(auto) ReactorDomain::begin() const
  {
    return containers.begin();
  }
  inline decltype(auto) ReactorDomain::end() const
  {
    return containers.end();
  }

  inline decltype(auto) ReactorDomain::begin()
  {
    return containers.begin();
  }
  inline decltype(auto) ReactorDomain::end()
  {
    return containers.end();
  }

} // namespace MC

#endif //__MC_REACTORDOMAIN_HPP__