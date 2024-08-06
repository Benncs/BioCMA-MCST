#ifndef __MC_PARTICLESHPP__
#define __MC_PARTICLESHPP__

#include "cmt_common/macro_constructor_assignment.hpp"
#include "common/execinfo.hpp"
#include "mc/prng/prng.hpp"
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
  template <typename T> class alignas(ExecInfo::cache_line_size) BaseParticle
  {
  public:
    BaseParticle() noexcept
        : current_container(0), current_domain(0), random_seed(0), id(0),
          status(CellStatus::IDLE), weight(0.){};

    explicit BaseParticle(double _weight) noexcept
        : current_container(0), current_domain(0), random_seed(0), id(0),
          status(CellStatus::IDLE), weight(_weight){};

    void clearState(MC::CellStatus _status = CellStatus::IDLE) noexcept
    {
      current_container = 0;
      current_domain = 0;
      random_seed = 0;
      id = 0;
      status = _status;
      weight = 0;



      clearData();
    }

    BaseParticle(const BaseParticle &p) = default; // Copy constructor
    BaseParticle &
    operator=(const BaseParticle &p) = default; // Copy assignment operator
    BaseParticle(BaseParticle &&p) noexcept = default;
    BaseParticle &
    operator=(BaseParticle &&p) noexcept = default; // Move assignment operator

    ~BaseParticle() = default;
    size_t current_container;
    size_t current_domain;
    size_t random_seed;
    uint32_t id;
    std::shared_ptr<MC::KPRNG> rng;
    MC::CellStatus status;
    double weight;
    T data;

    template <class Archive> void serialize(Archive &ar)
    {
      ar(current_container, current_domain, random_seed, id, status, weight);
    }

  private:
    void clearData() noexcept
    {

      if constexpr (requires { data.clear(); })
      {
        data.clear();
      }
      else if constexpr (requires { data.reset(); })
      {
        data.reset();
      }
      else
      {
        data.~T();
      }
    }
  };

  using Particles = BaseParticle<std::any>;

  
} // namespace MC

#endif