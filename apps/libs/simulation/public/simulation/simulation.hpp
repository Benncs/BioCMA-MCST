#ifndef __SIMULATIONS_UNIT_HPP__
#define __SIMULATIONS_UNIT_HPP__

#include <Kokkos_Core.hpp>
#include <Kokkos_ScatterView.hpp>
#include <biocma_cst_config.hpp>
#include <cassert>
#include <cma_read/reactorstate.hpp>
#include <common/common.hpp>
#include <common/kokkos_vector.hpp>
#include <cstddef>
#include <cstdint>
#include <mc/events.hpp>
#include <mc/particles/extra_process.hpp>
#include <mc/prng/prng.hpp>
#include <mc/unit.hpp>
#include <memory>
#include <optional>
#include <simulation/alias.hpp>
#include <simulation/feed_descriptor.hpp>
#include <simulation/probe.hpp>
#include <simulation/scalar_initializer.hpp>
#include <simulation/simulation_kernel.hpp>
// #include <simulation/kernel_update.hpp>
// #include <simulation/kernel_move.hpp>

// TODO Clean
static constexpr size_t trigger_const_particle_number = 1e6;

/**
 * @namespace Simulation
 * @brief Namespace that contains classes and structures related to simulation
 * handling
 */
namespace Simulation
{

  struct Dimensions
  {
    std::size_t n_species{};
    std::size_t n_compartment{};

    template <class Archive> void serialize(Archive& archive)
    {
      archive(n_species, n_compartment);
    }
  };

  class PreCalculatedHydroState;

  class ScalarSimulation;

  class SimulationUnit
  {
  public:
    SimulationUnit(std::unique_ptr<MC::MonteCarloUnit>&& _unit,
                   const ScalarInitializer& scalar_init,
                   std::optional<Feed::SimulationFeed> _feed = std::nullopt);

    ~SimulationUnit() = default;

    SimulationUnit(SimulationUnit&& other) noexcept;
    SimulationUnit(const SimulationUnit& other) = delete;
    SimulationUnit& operator=(SimulationUnit&& rhs) = delete;
    SimulationUnit& operator=(const SimulationUnit& rhs) = delete;

    std::unique_ptr<MC::MonteCarloUnit> mc_unit;

    [[nodiscard]] std::span<double> getCliqData() const;
    [[nodiscard]] std::optional<std::span<const double>> getCgasData() const;
    [[nodiscard]] Dimensions getDimensions() const noexcept;
    [[nodiscard]] std::span<const double> getContributionData() const;
    [[nodiscard]] const Simulation::Feed::SimulationFeed& get_feed() const;
    ProbeAutogeneratedBuffer& get_probes();

    void setVolumes(std::span<const double> volumesgas, std::span<const double> volumesliq) const;
    void setLiquidFlow(PreCalculatedHydroState* _flows_l);
    void setGasFlow(PreCalculatedHydroState* _flows_g);
    void set_probes(ProbeAutogeneratedBuffer&& _probes);

    // Simulation methods
    void cycleProcess(auto&& container, double d_t);

    void step(double d_t, const CmaRead::ReactorState& state) const;

    void reduceContribs(std::span<const double> data, size_t n_rank) const;

    [[deprecated("perf:not useful")]] void
    reduceContribs_per_rank(std::span<const double> data) const;

    void clearContribution() const noexcept;

    void update_feed(double t, double d_t, bool update_scalar = true) noexcept;

    [[nodiscard]] bool two_phase_flow() const;

    [[nodiscard]] std::size_t counter() const;

    [[nodiscard]] double& get_start_time_mut();

    [[nodiscard]] double& get_end_time_mut();

    template <class ListType, class CompartmentListType>
    inline void post_kernel_process(ListType& list,
                                    ListType& extra_list,
                                    CompartmentListType& local_compartments,
                                    size_t _waiting_allocation_particle);
    // Memory management
    void clear_mc();

    void reset();
    void post_init_compartments();

  private:
    // Attributes
    Kokkos::View<std::size_t, Kokkos::SharedSpace> internal_counter_dead;
    Kokkos::View<size_t, Kokkos::SharedSpace> waiting_allocation_particle;
    ProbeAutogeneratedBuffer probes;
    Simulation::Feed::SimulationFeed feed;

    PreCalculatedHydroState* flow_liquid; // TODO OPTI
    PreCalculatedHydroState* flow_gas;    // TODO OPTI
    bool const_number_simulation = true;
    bool is_two_phase_flow;
    double starting_time = 0.; // Not used within calculation, only for export purposes
    double end_time{};         // Not used within calculation, only for export purposes

    // Bounce methods to pimpl
    DiagonalView<ComputeSpace> get_kernel_diagonal();
    CumulativeProbabilityViewCompute get_kernel_cumulative_proba();
    kernelContribution get_kernel_contribution();
    [[nodiscard]] NeighborsView<ComputeSpace> get_kernel_neighbors() const;
    void set_kernel_contribs_to_host();

    void post_init_concentration(const ScalarInitializer& scalar_init);
    void post_init_concentration_functor(const ScalarInitializer& scalar_init);
    void post_init_concentration_file(const ScalarInitializer& scalar_init);

    struct pimpl_deleter
    {
      void operator()(ScalarSimulation*) const;
    };

    using pimp_ptr_t = std::unique_ptr<ScalarSimulation, pimpl_deleter>;
    pimp_ptr_t liquid_scalar;
    pimp_ptr_t gas_scalar;

    KernelInline::MoveInfo move_info;
  };

  inline ProbeAutogeneratedBuffer& SimulationUnit::get_probes()
  {
    return probes;
  }

  inline void SimulationUnit::setLiquidFlow(PreCalculatedHydroState* _flows_l)
  {
    flow_liquid = _flows_l;
  }

  inline void SimulationUnit::setGasFlow(PreCalculatedHydroState* _flows_g)
  {
    flow_gas = _flows_g;
  }

  void SimulationUnit::cycleProcess(auto&& container, double d_t)
  {

    PROFILE_SECTION("cycleProcess")
    auto& list = container.get_compute();
    auto& extra_list = container.extra_list;

    const size_t n_particle = list.size();


 
    auto& local_compartments = mc_unit->domain.data();
    auto& local_rng = mc_unit->rng;
    auto events = mc_unit->events;
    auto contribs = get_kernel_contribution();

    ContributionView contribs_scatter(contribs);
    Kokkos::deep_copy(waiting_allocation_particle, 0);

    this->move_info.cumulative_probability = get_kernel_cumulative_proba();
    this->move_info.diag_transition=get_kernel_diagonal();
    this->move_info.neighbors = mc_unit->domain.getNeighbors();

    Kokkos::parallel_for("mc_cycle_process",
                         Kokkos::RangePolicy<ComputeSpace>(0, n_particle),
                         KernelInline::Kernel(d_t,
                                              list,
                                              extra_list,
                                              internal_counter_dead,
                                              local_compartments,
                                              this->move_info,
                                              events,
                                              contribs_scatter,
                                              probes,
                                              waiting_allocation_particle,
                                              local_rng.random_pool));

    Kokkos::fence("fence_mc_cycle_process");

    Kokkos::Experimental::contribute(contribs, contribs_scatter);
    const std::size_t _waiting_allocation_particle = waiting_allocation_particle();
    post_kernel_process(list, extra_list, local_compartments, _waiting_allocation_particle);
  }

  template <class ListType, class CompartmentListType>
  inline void SimulationUnit::post_kernel_process(ListType& list,
                                                  ListType& extra_list,
                                                  CompartmentListType& local_compartments,
                                                  size_t _waiting_allocation_particle)
  {

    Kokkos::parallel_for(
        "update_compartment_number", extra_list.size(), KOKKOS_LAMBDA(const int i) {
          Kokkos::atomic_increment(
              &local_compartments(extra_list._owned_data(i).properties.current_container).n_cells);
        });

    static constexpr uint64_t minimum_dead_particle_removal = 100;
    const auto threshold =
        std::max(minimum_dead_particle_removal,
                 static_cast<uint64_t>(static_cast<double>(list.size()) *
                                       AutoGenerated::dead_particle_ratio_threshold));

    ;
    // Use commented line to remove dead when there is no particle.
    // This lead to stop exporting state for next time steps.
    //  if (_internal_counter > threshold || _internal_counter==list.size())
    if (internal_counter_dead() > threshold)
    {
#ifndef NDEBUG
      const auto old_size = list.size();
#endif
      list.remove_dead(internal_counter_dead());
#ifndef NDEBUG
      KOKKOS_ASSERT(list.size() == old_size - internal_counter_dead());
#endif
      internal_counter_dead() = 0;
    }

    list.insert(extra_list);
    extra_list.clear();
    const size_t n_new_alloc = _waiting_allocation_particle;
    if (n_new_alloc != 0)
    {
      const double new_weight = this->mc_unit->init_weight;
      list._spawn_alloc(n_new_alloc, new_weight);

      Kokkos::parallel_for(
          "add_new_alloc",
          Kokkos::RangePolicy<ComputeSpace>(0, n_new_alloc),
          KOKKOS_LAMBDA(const int) { Kokkos::atomic_increment(&local_compartments(0).n_cells); });
      Kokkos::fence("Fence cycle process ");
    }

    set_kernel_contribs_to_host();
  }

} // namespace Simulation

#endif //__SIMULATIONS_UNIT_HPP__
