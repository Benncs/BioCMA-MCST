// #include "mc/particles/mcparticles.hpp"
// #include "mc/particles/particles_list.hpp"
// #include <algorithm>
// #include <iostream>
// #include <mc/particles/particles_container.hpp>
// #include <omp.h>
// #include <stdexcept>
// #include <vector>

// namespace MC
// {
 
//   void ParticlesContainer::merge(Results& i_data)
//   {
//     auto &dead = i_data.index_in_dead_state;
//     auto &new_p = i_data.extra_process;

//     size_t new_size = new_p.size();
//     size_t dead_count = dead.size();

//     if (to_process.size() + (new_size - dead_count) > MAX_PARTICLE_BUFFER)
//     {
//       throw std::runtime_error("Overflow: particle list size exceeds limits");
//     }

//     auto new_p_iter = new_p.begin();
//     auto dead_iter = dead.begin();


//     for (;dead_iter<dead.end();++dead_iter)
//     {

//       if (new_p_iter != new_p.end())
//       {
//         to_process[*dead_iter] = std::move(*new_p_iter);
//         ++new_p_iter;
//       }
//       else
//       {
//         break;
//       }
//     }

//     dead.erase(dead_iter, dead.end());

//     to_process.insert(std::make_move_iterator(new_p_iter),
//         std::make_move_iterator(new_p.end()));
//     new_p.clear();
//   }

//   ParticlesContainer::ParticlesContainer(ParticlesContainer &&other) noexcept
//   {
//     if (this != &other)
//     {
//       to_process = std::move(other.to_process);
//     }
//   }

//   ParticlesContainer::ParticlesContainer(size_t capacity,
//                                          double weight) noexcept
//   {
//     to_process = ParticlesList(capacity, weight);
//     // init_extra(n_extra);
//   }

//   // void ParticlesContainer::init_extra(size_t n_extra)
//   // {
//   //   extras.resize(n_extra);
//   // }
// } // namespace MC