#ifndef __PARTICLES_CONTAINER_HPP__
#define __PARTICLES_CONTAINER_HPP__

#include <mc/particles/extra_process.hpp>
#include <mc/particles/particle_list.hpp>
#include <mc/particles/particle_model.hpp>
#include <mc/thread_private_data.hpp>
#include <variant>

namespace MC
{

  template <ParticleModel Model> class ParticlesContainer
  {

  public:
    using ComputeSpace = Kokkos::DefaultExecutionSpace::memory_space;
    using HostSpace = Kokkos::HostSpace::memory_space;
    DEFAULT_COPY_MOVE_AC(ParticlesContainer<Model>)

    ~ParticlesContainer() = default;

    auto &get_compute()
    {
      return to_process;
    }
    explicit ParticlesContainer(size_t capacity, double weight) noexcept
        : to_process(ParticlesList<ComputeSpace, Model>(capacity)),
          extra(capacity * 2)
    {

      to_process.init(weight);
    }

    auto &get_host()
    {
      ParticleList<ComputeSpace, Model>::migrate(to_process, host_process);
      return host_process;
    }

  private:
    Results<ComputeSpace, Model> extra;
    ParticleList<ComputeSpace, Model> to_process;
    ParticleList<HostSpace, Model> host_process;
  };

  using ContainerVariant = std::variant<ParticlesContainer<DefaultModel>>;

} // namespace MC

#endif