#ifndef __MC_PARTICLESHPP__
#define __MC_PARTICLESHPP__

#include "Kokkos_Assert.hpp"
#include "Kokkos_Macros.hpp"
#include "Kokkos_Printf.hpp"
#include "common/execinfo.hpp"
#include "common/has_serialize.hpp"
#include "common/kokkos_vector.hpp"
#include "models/model_monod.hpp"
#include <mc/particles/data_holder.hpp>
#include <mc/particles/particle_model.hpp>
#include <mc/prng/prng.hpp>
#include <type_traits>
#include <utility>

namespace MC
{

  template <ParticleModel _Model> class alignas(ExecInfo::cache_line_size) BaseParticle
  {
  public:
    using Model = _Model;

    KOKKOS_INLINE_FUNCTION explicit BaseParticle(double _weight = 0) noexcept
        : properties(_weight){};

    KOKKOS_INLINE_FUNCTION void clearState(MC::CellStatus _status = CellStatus::IDLE) noexcept
    {
      properties.reset();
      properties.status = _status;
    }

    KOKKOS_INLINE_FUNCTION void init(KPRNG globalrng)
    {
      data.init(properties, globalrng);
    }

    KOKKOS_INLINE_FUNCTION void
    update(double d_t,
           const LocalConcentrationView& concentration,
           MC::KPRNG globalrng)
    {
      properties.status = MC::CellStatus::IDLE;
      properties.hydraulic_time += d_t;
      properties.interdivision_time += d_t;
      data.update(d_t, properties, concentration, globalrng);
    }

    KOKKOS_INLINE_FUNCTION BaseParticle<_Model> division(KPRNG globalrng)
    {
      this->properties.status = CellStatus::IDLE;
      this->properties.interdivision_time = 0;
      auto prop_child = this->properties;
      prop_child.hydraulic_time = 0;
      auto p = data.division(prop_child, globalrng);

   

      return BaseParticle(std::move(prop_child), std::move(p));
    }

    KOKKOS_INLINE_FUNCTION void contribution(const ContributionView& contrib)
    {
      data.contribution(properties, contrib);
    }

    template <class Archive> void serialize(Archive& ar)
    {

      ar(properties);
      if constexpr (has_serialize<_Model, Archive>())
      {
        ar(data);
      }
    }

    alignas(ExecInfo::cache_line_size)  ParticleDataHolder properties;
    alignas(ExecInfo::cache_line_size)  _Model data{};

    KOKKOS_INLINE_FUNCTION BaseParticle(ParticleDataHolder&& props, _Model&& _model)
        : properties(props), data(std::move(_model))
    {
    }

    KOKKOS_INLINE_FUNCTION void fill_properties(
        SubViewtype subview,
        const Kokkos::Experimental::ScatterView<double**, Kokkos::LayoutRight, ComputeSpace>&
            spatial)
    {
      KOKKOS_ASSERT(subview.extent(0) == (2 + Model::get_number()));

      constexpr std::size_t size = Model::get_number();

      subview(size) = properties.hydraulic_time;
      subview(size + 1) = properties.interdivision_time;
      data.fill_properties(subview);
      {
        auto access_contribs = spatial.access();
        access_contribs(size, properties.current_container) += properties.hydraulic_time;
        access_contribs(size + 1, properties.current_container) += properties.interdivision_time;
        for (std::size_t i = 0; i < size; ++i)
        {
          access_contribs(i, properties.current_container) += subview(i);
        }
      }
    }
  };

  template <ParticleModel Model> using Particle = BaseParticle<Model>;

} // namespace MC

#endif
