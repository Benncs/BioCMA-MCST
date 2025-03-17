#include "mc/prng/prng.hpp"
#include "mc/traits.hpp"
#include <Kokkos_Core.hpp>
#include <Kokkos_ScatterView.hpp>
#include <common/kokkos_vector.hpp>
#include <cstdint>
#include <cstring>
#include <mc/unit.hpp>
#include <type_traits>
#include <utility>
#include <variant>

namespace
{
  /**
   * @brief Functor to count the number of particles in each container.
   *
   * This functor increments the count of particles in their respective containers.
   *
   * @tparam ListType The type of the particle list.
   * @tparam ScatterType The type of the scatter view for atomic updates.
   */
  template <typename ScatterType> struct NcellFunctor
  {
    ScatterType n_cells;             ///< Scatter view for counting particles in containers.
    MC::ParticlePositions positions; ///< List of particles.
    MC::ParticleStatus status;       ///< List of particles.

    /**
     * @brief Constructor for NcellFunctor.
     * @param _list The particle list.
     * @param _n_cells The scatter view for counting.
     */
    KOKKOS_FUNCTION
    NcellFunctor(MC::ParticlePositions _list, MC::ParticleStatus _status, ScatterType _n_cells)
        : n_cells(std::move(_n_cells)), positions(std::move(_list)), status(std::move(_status))
    {
    }

    /**
     * @brief Operator to update the particle container count.
     * @param i_particle The index of the particle being processed.
     */
    KOKKOS_FUNCTION void operator()(const int i_particle) const
    {
      auto access_ = n_cells.access();
      // access_(list._owned_data(i_particle).properties.current_container) += 1;
      if (status(i_particle) == MC::Status::Idle)
      {
        access_(positions(i_particle)) += 1;
      }
    }
  };

  /**
   * @brief Functor to initialize particle weights after creation.
   *
   * This functor updates the weight property of all particles in the list.
   *
   * @tparam ListType The type of the particle list.
   */
  struct PostInitFunctor
  {
    MC::ParticleWeigths weights; ///< List of particles.
    double new_weigth;        ///< The new weight assigned to each particle.

    /**
     * @brief Constructor for PostInitFunctor.
     * @param _list The particle list.
     * @param _new_weigth The new weight to be assigned to each particle.
     */
    PostInitFunctor(MC::ParticleWeigths _weights, double _new_weigth)
        : weights(std::move(_weights)), new_weigth(_new_weigth)
    {
    }

    /**
     * @brief Operator to update particle weights.
     * @param i_particle The index of the particle being processed.
     */
    KOKKOS_INLINE_FUNCTION void operator()(std::size_t i_particle) const
    {
      // list._owned_data(i_particle).properties.weight = new_weigth;
    }
  };

  template <ModelType Model> struct InitFunctor
  {
    explicit InitFunctor(MC::ParticlesContainer<Model> _list,
                         uint64_t min,
                         uint64_t max,
                         MC::KPRNG kprng)
        : particles(std::move(_list)), rng(std::move(kprng)), min_c(min), max_c(max) {

          };

    KOKKOS_INLINE_FUNCTION void operator()(const int i, double& local_mass) const
    {
      // // auto particle_rng = list.rng_instance;
      // auto p = MC::Particle<CurrentModel>(1.);
      // p.properties.id = i;
      // const uint64_t location = rng.uniform_u(min_c, max_c);
      particles.position(i) = rng.uniform_u(min_c, max_c);
      Model::init(rng.random_pool,i, particles.model);
      // p.properties.current_container = location;
      // p.init(rng);
      const double mass_i = Model::mass(i,particles.model);
      local_mass += mass_i;
      // list.set(i, std::move(p));
    }

    MC::ParticlesContainer<Model> particles;
    MC::KPRNG rng;
    uint64_t min_c;
    uint64_t max_c;
  };

} // namespace

namespace MC
{
  [[nodiscard]] uint64_t MonteCarloUnit::n_particle() const
  {
    return std::visit([](auto&& c) -> std::size_t { return c.n_particles(); }, container);
  }

  [[nodiscard]] std::vector<uint64_t> MonteCarloUnit::getRepartition() const
  {
    if (domain.getNumberCompartments() == 1)
    {
      return {n_particle()};
    }
    // Compute on-deman ncells as long as we don´t need it during simulation but only while
    // exporting. If ncells has to be known more often (id during iteration), n_cell view can be
    // allocated in mc_unit and incremented during kernel using atomic.
    Kokkos::View<uint64_t*, ComputeSpace> n_cells("n_cell", this->domain.getNumberCompartments());
    Kokkos::Experimental::ScatterView<uint64_t*> sn_cells(n_cells);
    std::visit(
        [&sn_cells](auto&& _container)
        {
          // auto& list = _container.get_compute();
          Kokkos::parallel_for("get_repartition",
                               _container.n_particles(),
                               NcellFunctor(_container.position, _container.status, sn_cells));
        },
        container);
    Kokkos::Experimental::contribute(n_cells, sn_cells);
    auto host = Kokkos::create_mirror_view_and_copy(HostSpace(), n_cells);
    std::vector<uint64_t> dist(n_cells.extent(0));
    std::memcpy(dist.data(), host.data(), host.size() * sizeof(uint64_t));

    return dist;
  }

  void post_init_weight(std::unique_ptr<MonteCarloUnit>& unit, double x0, double total_mass)
  {
    auto functor = [total_mass, x0, &unit](auto& container)
    {
      // TODO

      // const std::size_t n_p = container.n_particles();
      const double new_weight = (x0 * unit->domain.getTotalVolume()) / (total_mass);
      KOKKOS_ASSERT(new_weight > 0);
      using CurrentModel = typename std::remove_reference<decltype(container)>::type::UsedModel;
      if constexpr (ConstWeightModelType<CurrentModel>)
      {
        Kokkos::deep_copy(container.weights, new_weight);
      }
      else
      {
        static_assert(!ConstWeightModelType<CurrentModel>, "Multiple weights Not implemented yet");
        // Kokkos::parallel_for("mc_init_apply",
        //                    Kokkos::RangePolicy<ComputeSpace>(0, n_p),
        //                    PostInitFunctor(container.weights, new_weight));
      }
      // Kokkos::View<double, ComputeSpace> view_new_weight("view_new_weight");
      // Kokkos::deep_copy(view_new_weight, new_weight);
      
      unit->init_weight = new_weight;
    };

    std::visit(functor, unit->container);
  }

  void impl_init(double& total_mass,
                 uint64_t n_particles,
                 MonteCarloUnit& unit,
                 ContainerVariant&& container)
  {
    constexpr bool uniform_init = true;
    auto visitor = [&](auto&& container)
    {
      auto rng = unit.rng;

      // list.set_allocation_factor(AutoGenerated::default_particle_container_allocation_factor);
      const auto n_compartments = unit.domain.getNumberCompartments();
      uint64_t min_c = 0;
      uint64_t max_c = n_compartments;

      if (!uniform_init)
      {
        min_c = 0;
        max_c = 1;
      }

      /*
       * The mass of each cell in the reactor can be calculated after the model initialization.
       * To ensure that each cell has a unique weight based on the total mass, the following
       * formula is used: weight = XV / m_tot Where: XV  - represents a certain property or value
       * related to the cell (e.g., volume, particle count, etc.). m_tot - the total mass of the
       * cell or reactor, which needs to be determined first.
       *
       * In order to compute the weight correctly, the initialization process needs to be split
       * into two phases:
       * 1. The first phase is to calculate the total mass of the cell (m_tot).
       * 2. The second phase is to apply the newly calculated weight to the cell using the formula
       * above.
       *
       * This split ensures that the mass is determined before the weight, as the weight is
       * dependent on the total mass.
       */
      using CurrentModel = typename std::remove_reference<decltype(container)>::type::UsedModel;

      Kokkos::parallel_reduce("mc_init_first",
                              Kokkos::RangePolicy<Kokkos::DefaultExecutionSpace>(0, n_particles),
                              InitFunctor<CurrentModel>(container, min_c, max_c, rng),
                              total_mass);
      Kokkos::fence();
    };

    std::visit(visitor, container);
    unit.container = std::move(container);
  }

} // namespace MC