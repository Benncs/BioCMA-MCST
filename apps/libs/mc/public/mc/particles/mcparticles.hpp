#ifndef __MC_PARTICLESHPP__
#define __MC_PARTICLESHPP__

#include <Kokkos_Core.hpp>
#include <common/execinfo.hpp>
#include <common/has_serialize.hpp>
#include <mc/particles/data_holder.hpp>
#include <mc/particles/particle_model.hpp>
#include <mc/prng/prng.hpp>
#include <utility>

namespace MC
{

  /**
  @brief Represent Monte-Carlo particle, transport common and model specific propreties.

  The Model given as template  derived from ParticleModel therefore ensure compile time checking and
  slight more optimisation
   */
  template <ParticleModel _Model> class alignas(ExecInfo::cache_line_size) BaseParticle
  {
  public:
    using Model = _Model;

    /**
    @brief Default constructor

    _weight attribute can be given or not, if not given user has to ensure that it'll be set during
    main initialisation
     */
    KOKKOS_INLINE_FUNCTION explicit BaseParticle(double _weight = 0) noexcept
        : properties(_weight){};

    /**
    @brief Clear Particle State, if no state if given the newmade particle is in idle state
     */
    KOKKOS_INLINE_FUNCTION void clearState(MC::CellStatus _status = CellStatus::IDLE) noexcept
    {
      properties.reset();
      properties.status = _status;
    }

    /**
    @brief Bounce method to init particle model model
     */
    KOKKOS_INLINE_FUNCTION void init(KPRNG globalrng)
    {
      data.init(properties, globalrng);
    }

    /**
     * @brief Updates the particle state based on current concentrations.
     *
     * This method advances the particle's internal state.
     *
     * @param d_t Time step for the update.
     * @param concentration Local view of the current concentration.
     * @param globalrng Random number generator used for stochastic updates.
     */
    KOKKOS_INLINE_FUNCTION void
    update(double d_t, const LocalConcentrationView& concentration, MC::KPRNG globalrng)
    {
      properties.status = MC::CellStatus::IDLE;
      properties.hydraulic_time += d_t;
      properties.interdivision_time += d_t;
      data.update(d_t, properties, concentration, globalrng);
    }

    /**
     * @brief Performs particle division and returns a newly created particle.
     *
     * The properties of the new particle are determined by the selected model.
     * This method uses a random number generator to introduce stochasticity where applicable.
     *
     * @param globalrng Random number generator used during division.
     * @return BaseParticle<_Model> New particle with properties set according to the division
     * model.
     */
    KOKKOS_INLINE_FUNCTION BaseParticle<_Model> division(KPRNG globalrng)
    {
      this->properties.status = CellStatus::IDLE;
      this->properties.interdivision_time = 0;
      auto prop_child = this->properties;
      prop_child.hydraulic_time = 0;
      return BaseParticle(std::move(prop_child), std::move(data.division(prop_child, globalrng)));
    }

    /**
    @brief Bounce method to perform particle contribution
     */
    KOKKOS_INLINE_FUNCTION void contribution(const ContributionView& contrib)
    {
      data.contribution(properties, contrib);
    }

    /**
    @brief SerDe method
     */
    template <class Archive> void serialize(Archive& ar)
    {

      ar(properties);
      if constexpr (has_serialize<_Model, Archive>())
      {
        ar(data);
      }
    }

    /**
    @brief Constructor from existing model and properties

    This constructor ensure no copy by forwaring argument when possible
     */
    KOKKOS_INLINE_FUNCTION BaseParticle(ParticleDataHolder&& props, _Model&& _model)
        : properties(props), data(std::move(_model))
    {
    }

    /**
     * @brief Exports properties of the particle into arrays.
     *
     * This method allows a particle to write its properties into the provided arrays
     * to exportmodel-specific properties.
     * Records in both the local subview and the spatial ScatterView.
     *
     * @param subview Local subview where the particle's properties are stored.
     * @param spatial ScatterView for accumulating particle contributions across space.
     */
    KOKKOS_INLINE_FUNCTION void fill_properties(
        SubViewtype subview,
        const Kokkos::Experimental::ScatterView<double**, Kokkos::LayoutRight, ComputeSpace>&
            spatial)
    {
      KOKKOS_ASSERT(subview.extent(0) == (2 + Model::get_number()));

      const std::size_t size = Model::get_number();

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

    // Alignas ensure no falsesharing we accessing to properties only or data only
    alignas(ExecInfo::cache_line_size)
        ParticleDataHolder properties;                //< Particle's common properties
    alignas(ExecInfo::cache_line_size) _Model data{}; //< Particle's model
  };

  template <ParticleModel Model> using Particle = BaseParticle<Model>;

} // namespace MC

#endif
