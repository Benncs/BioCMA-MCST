#ifndef __MC_PARTICLESHPP__
#define __MC_PARTICLESHPP__

#include <any>
#include <cstddef>


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
    Particles()
        : current_container(0), current_domain(0), random_seed(0), id(0),
          status(CellStatus::IDLE), weight(0.){};
    explicit Particles(double _weight)
        : current_container(0), current_domain(0), random_seed(0), id(0),
          status(CellStatus::IDLE), weight(_weight){};

    void clearState(MC::CellStatus _status = CellStatus::IDLE)
    {
      current_container = 0;
      current_domain = 0;
      random_seed = 0;
      id = 0;
      status = _status;
      weight = 0;
      data.reset();
    }

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
    MC::CellStatus status;
    double weight;
    std::any data;
  };
} // namespace MC

#endif