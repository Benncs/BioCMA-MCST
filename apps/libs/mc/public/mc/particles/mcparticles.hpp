#ifndef __MC_PARTICLESHPP__
#define __MC_PARTICLESHPP__

#include "cmt_common/macro_constructor_assignment.hpp"
#include "common/execinfo.hpp"
#include <any>
#include <cstddef>
#include <cstdint>

namespace MC
{
  enum class CellStatus : char
  {
    IDLE,
    DEAD,
    CYTOKINESIS,
    OUT
  };

  class alignas(ExecInfo::cache_line_size) Particles
  {
  public:
    Particles() noexcept
        : current_container(0), current_domain(0), random_seed(0), id(0),
          status(CellStatus::IDLE), weight(0.){};

    explicit Particles(double _weight) noexcept
        : current_container(0), current_domain(0), random_seed(0), id(0),
          status(CellStatus::IDLE), weight(_weight){};

    void clearState(MC::CellStatus _status = CellStatus::IDLE) noexcept;

    Particles(const Particles &p) = default; // Copy constructor
    Particles &
    operator=(const Particles &p) = default; // Copy assignment operator
    Particles(Particles &&p) noexcept = default;
    Particles &
    operator=(Particles &&p) noexcept = default; // Move assignment operator

    ~Particles() = default;
    size_t current_container;
    size_t current_domain;
    size_t random_seed;
    uint32_t id;
    MC::CellStatus status;
    double weight;
    std::any data;

    template <class Archive> void serialize(Archive &ar)
    {
      ar(current_container, current_domain, random_seed, id, status, weight);
    }
  };

  inline void Particles::clearState(MC::CellStatus _status) noexcept
  {

    current_container = 0;
    current_domain = 0;
    random_seed = 0;
    id = 0;
    status = _status;
    weight = 0;
    data.reset();
  }
} // namespace MC

#endif