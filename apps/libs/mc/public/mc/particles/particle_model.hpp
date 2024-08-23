#ifndef __MC_PARTICLE_MODEL_HPP__
#define __MC_PARTICLE_MODEL_HPP__

#include <mc/particles/data_holder.hpp>

template <typename T>
concept ParticleModel = requires(T model, MC::ParticleDataHolder &p) {
  { model.init(p) } -> std::same_as<void>;
  { model.update(p) } -> std::same_as<void>;
  { model.division(p) } -> std::same_as<void>;
  { model.contribution(p) } -> std::same_as<void>;
};

struct DefaultModel
{
public:
  KOKKOS_FUNCTION void init(MC::ParticleDataHolder &p)
  {
    p.status = MC::CellStatus::IDLE;
  }

  KOKKOS_FUNCTION void update(MC::ParticleDataHolder &p)
  {
    p.current_container = (p.current_container + 1) % 5;
  }

  KOKKOS_FUNCTION void division(MC::ParticleDataHolder &p)
  {
    // Division logic
  }

  KOKKOS_FUNCTION void contribution(MC::ParticleDataHolder &p)
  {
    // Contribution logic
  }
};



#endif