#ifndef __MC_PARTICLESHPP__
#define __MC_PARTICLESHPP__

#include <any>
#include <cstddef>

#include <variant>

namespace MC
{
  enum class CellStatus
  {
    IDLE,
    DEAD,
    CYTOKINESIS
  };

  class Particles
  {
  public:
    Particles() = default;
    Particles(double _weight)
    {
      weight = _weight;
    };
    Particles(const Particles &p) = default; // Copy constructor
    Particles &
    operator=(const Particles &p) = default;     // Copy assignment operator
    Particles(Particles &&p) noexcept = default; // Move constructor
    Particles &
    operator=(Particles &&p) noexcept = default; // Move assignment operator

    size_t current_container;
    size_t current_domain;
    size_t random_seed;
    size_t id;
    MC::CellStatus status = CellStatus::IDLE;
    double weight;
    std::any data;
  };
} // namespace MC

#endif