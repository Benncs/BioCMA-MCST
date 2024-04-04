#ifndef __MC_PARTICLES_LIST_HPP__
#define __MC_PARTICLES_LIST_HPP__

#include <vector>

#include <mc/particles/mcparticles.hpp>

namespace MC
{
  class ParticlesList
  {
  public:
    ParticlesList() = default;

    ParticlesList(size_t capacity, double weight);

    ParticlesList(const ParticlesList &other) = delete;
    ParticlesList(ParticlesList &&other) noexcept;
    ~ParticlesList() = default;

    /*std::vector forward */
    void emplace_back(Particles &&p)
    {
      this->data.emplace_back(std::move(p));
    }

    // template <typename Iterator> void insert(Iterator &&begin, Iterator
    // &&end)
    // {
    //   this->data.insert(data.end(),begin, end);
    // }

    void insert(std::vector<MC::Particles> &&source)
    {
      data.insert(data.end(),
                  std::make_move_iterator(source.begin()),
                  std::make_move_iterator(source.end()));
    }


    
    inline void insert(MC::ParticlesList &&data)
    {
      // TODO
    }

    inline size_t size() const
    {
      return data.size();
    }

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

    auto &operator[](size_t i)
    {
      return data[i];
    }
    auto &operator[](size_t i) const
    {
      return data[i];
    }

    ParticlesList &operator=(const ParticlesList &other) = delete;

    ParticlesList &operator=(ParticlesList &&other) noexcept;

  private:
    std::vector<Particles> data;
  };
} // namespace MC

#endif