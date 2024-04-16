#include <mc/particles/particles_container.hpp>

namespace MC
{
  void ParticlesContainer::merge(size_t i_thread)
  {
    auto &dead = this->extras[i_thread].in_dead_state;
    auto &new_p = this->extras[i_thread].extra_process;
    int count = static_cast<int>(dead.size());
    int initial_size = count;
    if (count > 0)
    {
      for (auto &&i : new_p)
      {
        *dead[count - 1] = std::move(i);
        count--;
        if (count == 0)
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
                                         double weight) noexcept
  {
    to_process = ParticlesList(capacity, weight);
  }

  void ParticlesContainer::init_extra(size_t n_extra)
  {
    extras.resize(n_extra);
  }
} // namespace MC