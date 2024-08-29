#ifndef __PARTICLES_CONTAINER_HPP__
#define __PARTICLES_CONTAINER_HPP__

#include <mc/particles/extra_process.hpp>
#include <mc/particles/particle_list.hpp>
#include <mc/particles/particle_model.hpp>
#include <mc/thread_private_data.hpp>

namespace MC
{

  template <ParticleModel Model> class ParticlesContainer
  {

  public:
   
    using UsedModel = Model ;

    DEFAULT_COPY_MOVE_AC(ParticlesContainer)

    ~ParticlesContainer() = default;

    auto &get_compute()
    {
      return to_process;
    }
    explicit ParticlesContainer(size_t capacity) noexcept
        : to_process(capacity), host_process(), extra(capacity * 2)
    {

      // to_process.init(weight);
    }

    [[nodiscard]] inline size_t n_particle()const{return to_process.n_used_elements;}

    auto &get_host()
    {
      ParticleList<ComputeSpace, Model>::migrate(to_process, host_process);
      return host_process;
    }

    auto& get_extra(){return extra;}

    [[nodiscard]] size_t process_size() const
    {
      return to_process.size();
    }

  private:
    ParticleList<ComputeSpace, Model> to_process;
    ParticleList<HostSpace, Model> host_process;
    Results<ComputeSpace, Model> extra;
  };

 

} // namespace MC

#endif