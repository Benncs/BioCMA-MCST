#ifndef __MC_REACTORDOMAIN_HPP__
#define __MC_REACTORDOMAIN_HPP__

#include <cma_read/neighbors.hpp>
#include <common/common_types.hpp>
#include <mc/container_state.hpp>
#include <vector>

namespace MC
{
  class ReactorDomain
  {
  public:
    
    ReactorDomain() = default;
    ReactorDomain(ReactorDomain &&other) noexcept;
    ReactorDomain(std::span<double> volumes,
                  const Neighbors::Neighbors_const_view_t& _neighbors);


    ReactorDomain(const ReactorDomain &other) = delete;
    ~ReactorDomain() = default;

    ReactorDomain &operator=(const ReactorDomain &other) = delete;
    ReactorDomain &operator=(ReactorDomain &&other) noexcept;

    void setVolumes(std::span<double const> volumesgas,
                    std::span<double const> volumesliq);

    void setLiquidNeighbors(const Neighbors::Neighbors_view_t &data);

    // GETTERS
    auto &operator[](size_t i_c);
    auto &operator[](size_t i_c) const;

    [[nodiscard]] decltype(auto) begin() const;
    [[nodiscard]] decltype(auto) end() const;

    decltype(auto) begin();
    decltype(auto) end();

    [[nodiscard]] auto getNumberCompartments() const; 

    [[nodiscard]] double getTotalVolume() const;

    [[nodiscard]] const Neighbors::Neighbors_const_view_t  &getNeighbors() const;

    [[nodiscard]] std::span<const size_t> getNeighbors(size_t i) const;

    [[nodiscard]] std::vector<size_t> getDistribution() const;

    static ReactorDomain
    reduce(std::span<size_t> data, size_t original_size, size_t n_rank);

    std::span<ContainerState> data();

  private:
    double _total_volume=0.;
    size_t id=0;
    std::vector<ContainerState> containers;
    Neighbors::Neighbors_const_view_t neighbors;
  };

  inline std::span<ContainerState> ReactorDomain::data()
  {
    return containers;
  }

  inline std::span<const size_t> ReactorDomain::getNeighbors(size_t i) const
  {
    return neighbors.getRow(i);
  }

  inline const Neighbors::Neighbors_const_view_t &
  ReactorDomain::getNeighbors() const
  {
    
    return neighbors;
  }

  inline void ReactorDomain::setLiquidNeighbors(
      const Neighbors::Neighbors_view_t& data)
  {
 
    neighbors = data.to_const();
  }

  inline auto ReactorDomain::getNumberCompartments() const
  {
    return containers.size();
  }

  inline double ReactorDomain::getTotalVolume() const
  {
    return this->_total_volume;
  }

  inline auto &ReactorDomain::operator[](size_t i_c)
  {
    return this->containers[i_c];
  }
  inline auto &ReactorDomain::operator[](size_t i_c) const
  {
    return this->containers[i_c];
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