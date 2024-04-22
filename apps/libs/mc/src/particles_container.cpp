#include "mc/particles/particles_list.hpp"
#include <mc/particles/particles_container.hpp>
#include <stdexcept>

namespace MC
{
  void ParticlesContainer::merge(size_t i_thread)
  {
    auto &dead = this->extras[i_thread].in_dead_state;
    auto &new_p = this->extras[i_thread].extra_process;
    int dead_count = static_cast<int>(dead.size());
    if(to_process.size()+(new_p.size()-dead_count )>MAX_PARTICLE_BUFFER)
    {
      throw std::runtime_error("Overflow: particle list size exceeds limits") ;
      //TODO Try to catch this and load balance if MPI 
    }


    int initial_size = dead_count;
    
    if (dead_count > 0)
    {
      for (auto &&i : new_p)
      {
        *dead[dead_count - 1] = std::move(i);
        dead_count--;
        if (dead_count == 0)
        {
          break;
        }
      }
    }
    
    this->to_process.insert(new_p.begin() + initial_size, new_p.end());
    dead.clear();
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
                                         double weight,size_t n_extra) noexcept
  {
    to_process = ParticlesList(capacity, weight);
    init_extra(n_extra);
  }

  void ParticlesContainer::init_extra(size_t n_extra)
  {
    extras.resize(n_extra);
  }
} // namespace MC