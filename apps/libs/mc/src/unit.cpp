#include "mc/prng/prng.hpp"
#include <Kokkos_Core.hpp>
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
  template <typename ListType, typename ScatterType> struct NcellFunctor
  {
    ScatterType n_cells; ///< Scatter view for counting particles in containers.
    ListType list;       ///< List of particles.

    /**
     * @brief Constructor for NcellFunctor.
     * @param _list The particle list.
     * @param _n_cells The scatter view for counting.
     */
    KOKKOS_FUNCTION NcellFunctor(ListType _list, ScatterType _n_cells)
        : n_cells(std::move(_n_cells)), list(std::move(_list))
    {
    }

    /**
     * @brief Operator to update the particle container count.
     * @param i_particle The index of the particle being processed.
     */
    KOKKOS_FUNCTION void operator()(const int i_particle) const
    {
      auto access_ = n_cells.access();
      access_(list._owned_data(i_particle).properties.current_container) += 1;
    }
  };

  /**
   * @brief Functor to initialize particle weights after creation.
   *
   * This functor updates the weight property of all particles in the list.
   *
   * @tparam ListType The type of the particle list.
   */
  template <typename ListType> struct PostInitFunctor
  {
    ListType list;     ///< List of particles.
    double new_weigth; ///< The new weight assigned to each particle.

    /**
     * @brief Constructor for PostInitFunctor.
     * @param _list The particle list.
     * @param _new_weigth The new weight to be assigned to each particle.
     */
    PostInitFunctor(ListType _list, double _new_weigth)
        : list(std::move(_list)), new_weigth(_new_weigth)
    {
    }

    /**
     * @brief Operator to update particle weights.
     * @param i_particle The index of the particle being processed.
     */
    KOKKOS_INLINE_FUNCTION void operator()(std::size_t i_particle) const
    {
      list._owned_data(i_particle).properties.weight = new_weigth;
    }
  };

  template <typename ListType, typename CurrentModel> struct InitFunctor
  {
    explicit InitFunctor(ListType _list, uint64_t min, uint64_t max,MC::KPRNG kprng)
        : list(std::move(_list)), rng(std::move(kprng)),min_c(min), max_c(max) {

          };

    KOKKOS_INLINE_FUNCTION void operator()(const int i, double& local_mass) const
    {
      // auto particle_rng = list.rng_instance;
      auto p = MC::Particle<CurrentModel>(1.);
      p.properties.id = i;
      const uint64_t location = rng.uniform_u(min_c, max_c);
      p.properties.current_container = location;
      p.init(rng);
      const double mass_i = p.data.mass();
      local_mass += mass_i;
      list.set(i, std::move(p));
    }

    ListType list;
    MC::KPRNG rng;
    uint64_t min_c;
    uint64_t max_c;
  };

} // namespace

namespace MC
{
  [[nodiscard]] uint64_t MonteCarloUnit::n_particle() const
  {
    return std::visit([](auto&& c) -> std::size_t { return c.n_particle(); }, container);
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
          auto& list = _container.get_compute();
          Kokkos::parallel_for("get_repartition", list.size(), NcellFunctor(list, sn_cells));
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
      auto& list = container.get_compute();
      const std::size_t n_p = list.size();
      const double new_weight = (x0 * unit->domain.getTotalVolume()) / (total_mass);
      KOKKOS_ASSERT(new_weight > 0);
      Kokkos::View<double, ComputeSpace> view_new_weight("view_new_weight");
      Kokkos::deep_copy(view_new_weight, new_weight);
      Kokkos::parallel_for("mc_init_apply",
                           Kokkos::RangePolicy<ComputeSpace>(0, n_p),
                           PostInitFunctor(list, new_weight));
      unit->init_weight = new_weight;
    };

    std::visit(functor, unit->container);
  }

  void impl_init(double& total_mass,
                 uint64_t n_particles,
                 MonteCarloUnit& unit,
                 AutoGenerated::ContainerVariant&& container)
  {
    constexpr bool uniform_init = true;
    auto visitor = [&](auto&& c)
    {
      auto& list = c.get_compute();

      auto rng = unit.rng;
      auto particle_rng = list.rng_instance;

      list.set_allocation_factor(AutoGenerated::default_particle_container_allocation_factor);
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
       * To ensure that each cell has a unique weight based on the total mass, the following formula
       * is used: weight = XV / m_tot Where: XV  - represents a certain property or value related to
       * the cell (e.g., volume, particle count, etc.). m_tot - the total mass of the cell or
       * reactor, which needs to be determined first.
       *
       * In order to compute the weight correctly, the initialization process needs to be split into
       * two phases:
       * 1. The first phase is to calculate the total mass of the cell (m_tot).
       * 2. The second phase is to apply the newly calculated weight to the cell using the formula
       * above.
       *
       * This split ensures that the mass is determined before the weight, as the weight is
       * dependent on the total mass.
       */
      using CurrentModel = typename std::remove_reference<decltype(c)>::type::UsedModel;
      using ListType = typename std::remove_reference<decltype(list)>::type;
      // Kokkos::parallel_reduce(
      //     "mc_init_first",
      //     Kokkos::RangePolicy<Kokkos::DefaultExecutionSpace>(0, n_particles),
      //     KOKKOS_LAMBDA(const int i, double& local_mass) {
      //       auto p = Particle<CurrentModel>(1.);
      //       p.properties.id = i;
      //       const uint64_t location = rng.uniform_u(min_c, max_c);
      //       p.properties.current_container = location;
      //       p.init(particle_rng);
      //       const double mass_i = p.data.mass();
      //       local_mass += mass_i;
      //       list.set(i, std::move(p));
      //     },
      //     total_mass);

       Kokkos::parallel_reduce(
          "mc_init_first",
          Kokkos::RangePolicy<Kokkos::DefaultExecutionSpace>(0, n_particles),
          InitFunctor<ListType, CurrentModel>(list,min_c,max_c,rng),
          total_mass);
      Kokkos::fence();
    };

    std::visit(visitor, container);
    unit.container = std::move(container);
  }

} // namespace MC