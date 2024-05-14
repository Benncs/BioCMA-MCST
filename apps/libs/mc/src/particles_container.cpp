#include "mc/particles/particles_list.hpp"
#include <algorithm>
#include <mc/particles/particles_container.hpp>
#include <stdexcept>
// #include <iostream>
namespace MC
{
  void ParticlesContainer::merge(ThreadPrivateData &i_data)
  {
    auto &dead = i_data.in_dead_state;
    auto &new_p = i_data.extra_process;
    int new_size = static_cast<int>(new_p.size());
    int dead_count = static_cast<int>(dead.size());
    if (to_process.size() + (new_size - dead_count) > MAX_PARTICLE_BUFFER)
    {
      throw std::runtime_error("Overflow: particle list size exceeds limits");
      // TODO Try to catch this and load balance if MPI
    }

    int initial_size = dead_count;

    if (dead_count > 0)
    {
      for (auto &&i : new_p)
      {
        *dead[dead_count - 1] = std::move(i);
        dead[dead_count - 1] = nullptr;
        dead_count--;

        if (dead_count == 0)
        {
          break;
        }
      }

      dead.erase(std::remove_if(dead.begin(),
                                dead.end(),
                                [](auto *ptr) { return ptr == nullptr; }),
                 dead.end());
    }
    if (new_p.begin() + initial_size < new_p.end())
    {
      this->to_process.insert(new_p.begin() + initial_size, new_p.end());
    }

    new_p.clear();
  }

  ParticlesContainer::ParticlesContainer(ParticlesContainer &&other) noexcept
  {
    if (this != &other)
    {
      to_process = std::move(other.to_process);
    }
  }

  ParticlesContainer::ParticlesContainer(size_t capacity,
                                         double weight) noexcept
  {
    to_process = ParticlesList(capacity, weight);
    // init_extra(n_extra);
  }

  // void ParticlesContainer::init_extra(size_t n_extra)
  // {
  //   extras.resize(n_extra);
  // }
} // namespace MC