#ifndef __MC_PARTICLES_LIST_HPP__
#define __MC_PARTICLES_LIST_HPP__

#include <vector>

#include <mc/particles/mcparticles.hpp>

namespace MC
{
  class ParticlesList
  {
  public:
    explicit ParticlesList() = default;

    explicit ParticlesList(size_t capacity, double weight);

    explicit ParticlesList(const ParticlesList &other) = delete;
    explicit ParticlesList(ParticlesList &&other) noexcept;
    ~ParticlesList() = default;

    /*std::vector forward */
    void emplace_back(Particles &&p);

    void insert(std::vector<MC::Particles> &&source) noexcept;

    template <typename IT> void insert(IT &&begin, IT &&end) noexcept
    {
      data.insert(data.end(),
                  std::make_move_iterator(begin),
                  std::make_move_iterator(end));
    }

    size_t size() const noexcept;

    decltype(auto) begin() const
    {
      return data.begin();
    }
    decltype(auto) end() const
    {
      return data.end();
    }

    decltype(auto) begin()
    {
      return data.begin();
    }
    decltype(auto) end()
    {
      return data.end();
    }

    auto &operator[](size_t i);

    auto &operator[](size_t i) const;

    ParticlesList &operator=(const ParticlesList &other) = delete;

    ParticlesList &operator=(ParticlesList &&other) noexcept;

  private:
    std::vector<Particles> data;
  };

  inline void ParticlesList::emplace_back(Particles &&p)
  {
    this->data.emplace_back(std::move(p));
  }

  inline size_t ParticlesList::size() const noexcept
  {
    return data.size();
  }

  inline auto &ParticlesList::operator[](size_t i)
  {
    return data[i];
  }

  inline auto &ParticlesList::operator[](size_t i) const
  {
    return data[i];
  }

  inline void
  ParticlesList::insert(std::vector<MC::Particles> &&source) noexcept
  {
    data.insert(data.end(),
                std::make_move_iterator(source.begin()),
                std::make_move_iterator(source.end()));
  }
} // namespace MC

#endif